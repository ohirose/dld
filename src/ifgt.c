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
#include<math.h>
#include<assert.h>
#include<time.h>
#include"kmeans.h"
#include"kdtree.h"
#include"util.h"
#include"ifgt.h"

#define SQ(x) ((x)*(x))
#define TREE_MAXDEPTH 256

/* (state 1) push: S[++n]=p   <==> n++;    S[n]=p; */
/* (state 2) push: S[++n]=p++ <==> n++;p++;S[n]=p; */
/* (state *) pop:  S[n--]=p   <==> p=S[n]; n--;    */
void mchoose(int *A, int *M, int *S, int N, int K){
  int k,p=1,n=-1,m=0,j;
  /* initialization */
  for(k=0;k<K;k++) A[k]=0;
  S[++n]=p;

  if(N==0){*M=1;return;}

  tag:
  /* state 1 */
  while(n<N-1) S[++n]=p;
  for(j=0;j<N;j++) A[S[j]-1+m*K]++;
  m++; p=S[n--];

  /* state 2 */
  while(p==K&&n>=0) p=S[n--];
  if   (p==K&&n< 0) {*M=m;return;}
  if   (p<K)        {S[++n]=++p; goto tag;}
}

int nCk(int n, int k){
  int d,r=1; if(k>n) return 0;
  for(d=1;d<=k;d++){r*=n--;r/=d;}
  return r;
}


static double bstar(double a, int p, double h){return (a+sqrt(SQ(a)+2*p*SQ(h)))/2.0;}

static double delta(int p, double a, double b, double h){
  int q; double h2=SQ(h),c=2*a*b/h2; double val=1.0; 
  assert(p>=0); assert(a>=0&&b>=0);
  for(q=1;q<=p;q++) val*=c/(double)q;
  return val*exp(-SQ(a-b)/h2);
}

static double range(const double *X, int D, int N){
  double min=1e250,max=-1e250; int i,DN=D*N;
  for(i=0;i<DN;i++) min=fmin(min,X[i]);
  for(i=0;i<DN;i++) max=fmax(max,X[i]);
  return sqrt(D)*(max-min);
}

typedef enum ifgtmode{
  IFGT_MODE_TAYLOR    =0,
  IFGT_MODE_NEIGHBOR  =1,
  IFGT_MODE_DIRECT    =2,
} ifgtmode;

static void setting_init(int *mode /*O*/,int *K/*IO*/,int *Q/*IO*/,double *r/*O*/, int N, int D, double h, double e, double R){
  int k,p,P=*Q,KLIM=*K; int flag=0; double rx,len,n,b,c,val,enn,cmin=1e250;

  len =h/(R/sqrt(D));
  enn =pow(4*len,D)*N;
  val =e<1.0?h*sqrt(-log(e)):1e250;
  *r  =fmin(R,val);

  for(k=1;k<=KLIM;k++){rx=pow(k,-1.0/D); n=k; //n=fmin(k,pow(*r/rx,D));
    for(p=1;p<P;p++){
      b=fmin(bstar(rx,p,h),*r+rx);
      if(delta(p,rx,b,h)<e) break;
    }
    c=log(k)+(1+n)*nCk(p-1+D,D);
    if(p!=P&&c<cmin){cmin=c;*K=k;*Q=p;flag=1;}
  } 
  *K=flag?*K:0; 

  if   (enn>=cmin)  *mode=N>4.0*cmin?IFGT_MODE_TAYLOR  :IFGT_MODE_DIRECT;
  else/*enn< cmin*/ *mode=N>4.0*enn ?IFGT_MODE_NEIGHBOR:IFGT_MODE_DIRECT;

  return;
}

static void setting_update(int *mode /*O*/, int *Q/*IO*/, int K,  double rx, int N, int D, double h, double e, double R){
  int p,P=*Q; double enn,r,n,b,c,val;

  val =e<1.0?h*sqrt(-log(e)):1e250;
  enn =pow(4*h*sqrt(D)/R,D)*N;
  r   =fmin(R,val);
  n   =K; //fmin(K,pow(r/rx,D));

  for(p=1;p<P;p++){b=fmin(bstar(rx,p,h),r+rx); if(delta(p,rx,b,h)<e) break;}
  c=log(K)+(1+n)*nCk(p-1+D,D); *Q=p;
  
  if     (N>4.0*c&&*Q!=P) *mode=IFGT_MODE_TAYLOR;
  else if(N>4.0*enn)      *mode=IFGT_MODE_NEIGHBOR;
  else                    *mode=IFGT_MODE_DIRECT;

  return;
}

inline static double dist(const double *x, const double *y, int D){
  int d; double val=0;
  for(d=0;d<D;d++) val+=SQ(x[d]-y[d]);
  return sqrt(val);
}

inline static double dist2(const double *x, const double *y, int D){
  int d; double val=0;
  for(d=0;d<D;d++) val+=SQ(x[d]-y[d]);
  return val;
}

inline static double cterm(const double *z, const double *c, const int *a, int D, double h){
  int d; double val=1.0;
  for(d=0;d<D;d++) val*=pow((z[d]-c[d])/h,a[d]);
  return val;
}

/* make K and P large enough */
ifgtwork * ifgtwork_alloc(int M, int N, int D, int K, int P){
  ifgtwork *w; int sz=M>N?M:N,nt=nCk(P+D,D); int sd=sizeof(double),si=sizeof(int);

  w=calloc(1,sizeof(ifgtwork)); w->buf=calloc(2*sz,sd);
  w->L =calloc(P+1, si); w->U =calloc(D*K,   sd); w->g1=calloc(6*sz,si);
  w->A =calloc(D*nt,si); w->rx=calloc(1*K,   sd); w->g2=calloc(1*sz,si);
  w->B =calloc(1*nt,sd); w->G =calloc(1*sz,  si); w->g3=calloc(1*sz,si);
  w->C =calloc(K*nt,sd); w->T =calloc(3*sz+1,si); w->S =calloc2i(sz,TREE_MAXDEPTH);
  w->t =calloc(1*sz,si); w->a =calloc(1*sz,  si); 

  return w;
} 
  
void ifgt(
       double         *f,     /*  O  |   M   | resulting values              */
       ifgtwork       *w,     /* I/W |       | working memory                */
       const double   *Y,     /*  I  | D x M | input matrix Y                */
       const double   *X,     /*  I  | D x N | input matrix X                */ 
       const double   *q,     /*  I  |   N   | weight vector q               */
       int             D,     /*  I  | const.| dimension                     */
       int             M,     /*  I  | const.| #points in Y                  */
       int             N,     /*  I  | const.| #points in X                  */ 
       int             K,     /*  I  | const.| upper bound of #clusters      */
       int             P,     /*  I  | const.| upper bound of #trancations   */
       double          h,     /*  I  | const.| band width of Gauss func.     */
       double          e      /*  I  | const.| error bound                   */
     ){

  int d,j,k,l,n,m,p,Q=P,size,mode; double val,r,rad,dnk,dmk,snk,smk,h2=SQ(h),rxm=0;
  int Dl,Dn,Dk; double *rx,*U,*B,*C,*dbuf; int *t,*a,**S,*L,*A,*G,*tree,*nums,*ibuf,*stack;
  int nlp=10;
  
  /* working memory */
  L=w->L; A=w->A; B=w->B; C=w->C; U=w->U; rx=w->rx; G=w->G; tree=w->T; S=w->S; 
  t=w->t; a=w->a; dbuf=w->buf; ibuf=w->g1; stack=w->g2; nums=w->g3;

  /* initial parameters */
  val=e<1.0?sqrt(-log(e)):1e250; r=fmin(range(X,D,N),h*val);

  /* parameter estimation */
  setting_init(&mode,&K,&P,&r,N,D,h,e,range(X,D,N));

  switch(mode){
    case IFGT_MODE_TAYLOR:   goto taylor;   break;
    case IFGT_MODE_NEIGHBOR: goto neighbor; break;
    case IFGT_MODE_DIRECT:   goto direct;   break;
    default: printf("Error in IFGT mode. Abort.\n"); exit(EXIT_FAILURE); 
  }

  taylor:
  /* clustering */
  kmeans(G,U,nums,X,N,K,D,nlp);
  for(k=0;k<K;k++) rx[k]=0;
  for(n=0;n<N;n++){k=G[n];rx[k]=fmax(rx[k],dist(X+D*n,U+D*k,D));}
  for(k=0;k< K;k++) rxm=k==0?rx[k]:fmax(rx[k],rxm);

  setting_update(&mode,&Q,K,rxm,N,D,h,e,range(X,D,N)); P=Q;

  /* main computation for X */
  /* preparation */
  for(p=0;p<=P;p++) L[p]=nCk(p-1+D,D);
  for(p=0;p< P;p++) mchoose(A+D*L[p],&size,stack,p,D);
  for(l=0;l<L[P];l++){B[l]=1.0;Dl=D*l;for(d=0;d<D;d++)for(j=0;j<A[d+Dl];j++)B[l]*=2.0/(j+1);}
  for(l=0;l<L[P];l++)for(k=0;k<K;k++) C[k+K*l]=0;

  /* computation of C */
  for(n=0;n<N;n++){k=G[n];Dn=D*n;Dk=D*k;dnk=dist(X+Dn,U+Dk,D);snk=q[n]*exp(-SQ(dnk)/h2);rad=r+rx[k];
    for(p=0;p<P;p++){
      for(l=L[p];l<L[p+1];l++) C[k+K*l]+=B[l]*snk*cterm(X+Dn,U+Dk,A+D*l,D,h);
      if(delta(p,dnk,fmin(bstar(dnk,p,h),rad),h)<e) break;
    } 
  } 
  /* main computation for Y */
  #ifndef OMP_FOR_CMAES
    #pragma omp parallel for private(k) private(dmk) private(smk) private(p) private(l)
  #endif
  for(m=0;m<M;m++){f[m]=0;
    for(k=0;k<K;k++){dmk=dist(Y+D*m,U+D*k,D);if(dmk>r+rx[k]) continue; smk=exp(-SQ(dmk)/h2);
      for(p=0;p<P;p++)for(l=L[p];l<L[p+1];l++) f[m]+=C[k+K*l]*smk*cterm(Y+D*m,U+D*k,A+D*l,D,h);
    } 
  }
  return; 

  neighbor:
  kdtree(tree,ibuf,dbuf,X,D,N); rad=4.0*h; 
  #ifndef OMP_FOR_CMAES
    #pragma omp parallel for 
  #endif
  for(m=0;m<M;m++){a[m]=t[m]=f[m]=0;
    do{
      eballsearch_next(&a[m],S[m],&t[m],Y+D*m,rad,X,tree,D,N);
      if(a[m]>=0) f[m]+=q[a[m]]*exp(-dist2(X+D*a[m],Y+D*m,D)/h2);
    } while(t[m]);
  }
  return;

  direct:
  #ifndef OMP_FOR_CMAES
    #pragma omp parallel for private (n)
  #endif
  for(m=0;m<M;m++){f[m]=0;for(n=0;n<N;n++) f[m]+=q[n]*exp(-dist2(X+D*n,Y+D*m,D)/h2);}

}

