// Copyright (c) 2017 Osamu Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"util.h"
#include"cmaes.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int     dsyev_         (char *jobz,char *uplo,int *n,double *a,int *lda,double *w,double *work,int *lwork,int *info);
void    init_genrand   (unsigned long s);
double  genrand_res53  (void);

double norm(double *v, int D){int d;double val=0; for(d=0;d<D;d++) val+=SQ(v[d]); return sqrt(val);}

double variance (double *x, int N){
  int n; double m=0,v=0; assert(N>0);
  for(n=0;n<N;n++){m+=x[n];}       m/=N;
  for(n=0;n<N;n++){v+=SQ(x[n]-m);} v/=N-1;
  return v;
}

int normal(double *v, int nv){
  int i,n=nv/2,odd=(nv%2==1)?1:0; double u1,u2;
  if(n)for(i=0;i<n;i++){u1=genrand_res53();u2=genrand_res53();
    v[2*i  ]=sqrt(-2*log(u1))*cos(2*M_PI*u2);
    v[2*i+1]=sqrt(-2*log(u1))*sin(2*M_PI*u2);
  } 
  if(odd){u1=genrand_res53();u2=genrand_res53();v[nv-1]=sqrt(-2*log(u1))*cos(2*M_PI*u2);}
  return 0;
}

void print_cma_stats(int epoch, double score){
  int i;
  if(epoch)for(i=0;i<6;i++) printf("\033[F\033J");
  printf("\n  CMA-ES optimization:\n");
  printf("\n    Epoch: %d",epoch+1);
  printf("\n    Score: %lf\n\n", score);
  return;
}

void cma_size_init(cma_size *sz, int D){
  sz->D=D; 
  sz->N=4+floor(3*log(D));
}

cma_prms * cma_prms_alloc(cma_size sz){
  cma_prms *p=calloc(1,sizeof(cma_prms));
  p->w=calloc(sz.N,sizeof(double));
  return p;
}

void cma_prms_init(cma_prms *p, cma_size s){
  int n,N,D,mu; double u=0,v=0,u1=0,u2=0,v1=0,v2=0,*w=p->w;
  double csgm,dsgm,cc,c1,cmu,amu,amue,apd,mue,muen,acov=2.0;

  D=s.D; N=s.N; mu=floor(N/2);
  for(n=0;n<N;n++) w[n]=log(0.5*(1+N))-log(1+n);
  for(n=0;n<N;n++)if(n<mu) u1+=w[n];      else u2+=w[n];
  for(n=0;n<N;n++)if(n<mu) v1+=w[n]*w[n]; else v2+=w[n]*w[n];

  mue  =SQ(u1)/v1; u=sqrt((mue-1)/(D+1))-1;
  muen =SQ(u2)/v2; v=acov*(mue-2+1/mue)/(SQ(D+2)+0.5*acov*mue);

  assert(u1>=0); assert(mue>=1.0);  assert(muen>=0);
  assert(u2 <0); assert(mue<=mu);

  cc  =(4+mue/D)/(D+4+2*mue/D); csgm =(mue+2)/(D+mue+5); 
  c1  =acov/(SQ(D+1.3)+mue);    dsgm =1+2*MAX(0,u)+csgm; 
  cmu =MIN(v,1-c1);

  assert(csgm>0); assert(cc>0); assert(c1>0); assert(cmu>=0);
  assert(csgm<1); assert(cc<1); assert(c1<1); assert(cmu< 1);
  assert(dsgm>1);

  amu  =1+c1/cmu;            assert(amu >=0);
  amue =1+2*muen/(2+mue);    assert(amue>=0);
  apd  =(1-c1-cmu)/(D*cmu);  assert(apd >=0);

  for(n=0;n<N;n++)w[n]=(w[n]>=0)?(w[n]/u1):MIN(apd,MIN(amu,amue))*(w[n]/(-u2));

  p->csgm=csgm; p->c1=c1; p->cmu=cmu; p->mue=mue;
  p->dsgm=dsgm; p->cc=cc; p->mu=mu;
}

cma_stat *cma_stat_alloc(cma_size sz){
  int D=sz.D,N=sz.N,sod=sizeof(double),sos=sizeof(sortbox);
  cma_stat *s=calloc(1,sizeof(cma_stat));
  s->f   =calloc(N,sod); s->x =calloc(D*N,  sod);
  s->wc  =calloc(N,sod); s->y =calloc(D*N,  sod);
  s->L   =calloc(D,sod); s->z =calloc(D*N,  sod);
  s->yw  =calloc(D,sod); s->C =calloc(D*D,  sod);
  s->m   =calloc(D,sod); s->Cr=calloc(D*D,  sod);
  s->q   =calloc(D,sod); s->V =calloc(D*D,  sod);
  s->psgm=calloc(D,sod); s->A =calloc(D*D,  sod);
  s->pc  =calloc(D,sod); s->b =calloc(D*D+1,sod);
  s->sb  =calloc(N,sos);
  return s;
}

void cma_stat_init(cma_stat *s, cma_size sz, const double *m, const double sgm){
  int i,j,d,D=sz.D;
  s->sgm=sgm; 
  for(i=0;i<D;i++)for(j=0;j<D;j++) s->C[i+D*j]=i==j?1:0;
  for(d=0;d<D;d++) s->psgm[d]=0;
  for(d=0;d<D;d++) s->pc  [d]=0;
  for(d=0;d<D;d++) s->m   [d]=m[d];
}

void cma_optimize(
  double          *  xopt,                                /*  O  |  D  | x which optimizes f(x)              */
  cma_stat        *  stat,                                /*  W  |  *  | state: working memory for cmaes     */
  void            ** work,                                /*  W  |  N  | working memory used inside 'fptr'   */
  const void      ** args,                                /*  I  |  N  | arguments of f(x) except x and size */
  double          (*fptr)(double*,int,void*,const void*), /*  I  |  1  | f(x); x,size,work,arg               */
  const cma_prms  *  prms,                                /*  I  |  1  | parameters                          */
  const cma_conf     conf,                                /*  I  |  1  | configuration                       */
  const cma_size     size,                                /*  I  |  1  | D, N                                */
  const int          mode                                 /*  I  |  1  | maximize (1) or minimize (0)        */
  ){
  /* NOTE: 'args' and 'work' are defined as 1D arrays for pallarellization */
  int i,j,k,l,d,n,g,D,N,info,nlp,lb; double *w,*f,*m,*C,*x,*y,*z,*yw,*wc,*q,*psgm,*pc,*Cr,*L,*A,*V,*b;
  double e,a1,a2,csgm,dsgm,hsgm,dlt,c1,cc,cmu,mu,mue,sumw=0; char jobz='V',uplo='L'; sortbox *sb;
  double diff,tol,curr=0,prev=0;

  N =size.N; csgm=prms->csgm; c1=prms->c1; cmu=prms->cmu; mu=prms->mu;
  D =size.D; dsgm=prms->dsgm; cc=prms->cc; mue=prms->mue; w =prms->w;
  lb=D*D+1; e=sqrt(D)*(1-1/(4*D)+1/(21*D*D)); nlp=conf.epoc; tol=conf.tol;

  f=stat->f; x=stat->x; wc=stat->wc; V =stat->V; psgm=stat->psgm; 
  m=stat->m; y=stat->y; yw=stat->yw; A =stat->A; 
  L=stat->L; z=stat->z; pc=stat->pc; b =stat->b;
  q=stat->q; C=stat->C; Cr=stat->Cr; sb=stat->sb;

  /* initialization */
  init_genrand(conf.seed);
  for(n=0;n<N;n++) sumw+=w[n];

  for(g=0;g<nlp;g++){
    /* sample generation */
    normal(z,N*D);
    for(i=0;i<D;i++)for(j=0;j<D;j++) V[i+D*j]=C[i+D*j];
    dsyev_(&jobz,&uplo,&D,V,&D,L,b,&lb,&info); assert(!info);
    for(d=0;d<D;d++) L[d]=sqrt(L[d]);
    for(i=0;i<D;i++)for(j=0;j<D;j++) A[i+D*j]=V[i+D*j]*L[j];
    for(n=0;n<N;n++)for(d=0;d<D;d++){y[d+D*n]=0;for(k=0;k<D;k++) y[d+D*n]+=A[d+D*k]*z[k+D*n];}
    for(n=0;n<N;n++)for(d=0;d<D;d++) x[d+D*n]=m[d]+stat->sgm*y[d+D*n];

    /* function evaluation */
    #ifdef OMP_FOR_CMAES
      #pragma omp parallel for
    #endif
    for(n=0;n<N;n++) f[n]=fptr(x+D*n,D,(work==NULL?NULL:work[n]),(args==NULL?NULL:args[n]));

    /* selection and recombination */
    prepare_sortbox(sb,f,N); 
    if(mode) qsort(sb,N,sizeof(sortbox),cmp_sortbox_d);
    else     qsort(sb,N,sizeof(sortbox),cmp_sortbox_a);
    for(d=0;d<D;d++) yw[d]=0;
    for(d=0;d<D;d++)for(n=0;n<mu;n++){l=sb[n].idx; yw[d]+=w[n]*y[d+D*l];}
    for(d=0;d<D;d++) m[d]+=stat->sgm*yw[d];

    /* step-size control */
    for(i=0;i<D;i++)for(j=0;j<D;j++) Cr[i+D*j]=0;
    for(i=0;i<D;i++)for(j=0;j<D;j++)for(k=0;k<D;k++) Cr[i+D*j]+=V[i+D*k]*V[j+D*k]/L[k];
    for(d=0;d<D;d++){q[d]=0;for(i=0;i<D;i++) q[d]+=Cr[d+D*i]*yw[i];}
    for(d=0;d<D;d++) psgm[d]=(1-csgm)*psgm[d]+sqrt(csgm*(2-csgm)*mue)*q[d];
    for(d=0;d<D;d++) stat->sgm*=exp((csgm/dsgm)*(-1+(norm(psgm,D)/e)));

    /* covariance matrix adaptation */
    /* pc */
    a1=norm(psgm,D)/sqrt(1-pow(1-csgm,2*(g+1))); a2=e*(1.4+2/(D+1));
    hsgm=a1<a2?1.0:0.0; dlt=(1-hsgm)*cc*(2-cc); assert(hsgm<=1.0);
    for(d=0;d<D;d++) pc[d]=(1-cc)*pc[d]+hsgm*sqrt(cc*(2-cc)*mue)*yw[d];
    /* wc */
    for(n=0;n<N;n++){l=sb[n].idx;
      for(d=0;d<D;d++) q[d]=0;
      for(d=0;d<D;d++)for(k=0;k<D;k++) q[d]+=Cr[d+D*k]*y[k+D*l];
      //wc[n]=w[n]*(w[n]>=0?1.0:D/norm(q,D)); /* negative weights are not supproted */
      wc[n]=w[n]*(w[n]>=0?1.0:0.0); 
    }
    /* C  */
    for(i=0;i<D;i++)for(j=0;j<D;j++) C[i+D*j]*=1.0+c1*dlt-c1-cmu*sumw;
    for(i=0;i<D;i++)for(j=0;j<D;j++) C[i+D*j]+=c1*pc[i]*pc[j];
    for(n=0;n<N;n++){l=sb[n].idx;for(i=0;i<D;i++)for(j=0;j<D;j++) C[i+D*j]+=cmu*wc[n]*y[i+D*l]*y[j+D*l];}

    /* check convergence */
    prev=curr; curr=variance(f,N); diff=fabs(curr-prev);
    if(conf.verb) print_cma_stats(g,f[sb[0].idx]);
    if(g>1&&diff<tol) break;
  } 
  /* output */
  for(d=0;d<D;d++) xopt[d]=x[d+D*sb[0].idx];

  return;
}

