// Copyright (c) 2016-2017 Osamu Hirose
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
#include<string.h>
#include<assert.h>
#include<ctype.h>
#include<math.h>
#include"util.h"
#include"mat.h"
#include"error.h"
#include"ifgt.h"
#include"nystrom.h"
#include"dld.h"

/* macros for ifgt */
#define KLIM 1000
#define PLIM 20

enum dldsflags{
  FLAG_DEFORM         = (1 << 0),
  FLAG_SIMILAR        = (1 << 1),
  FLAG_IFGT           = (1 << 2),
  FLAG_NYSTROM        = (1 << 3),
  FLAG_LOCAL          = (1 << 4),
  FLAG_ADAPT          = (1 << 5),
  FLAG_GEMAN          = (1 << 6)
};

enum dldstate{
  STATE_START         = (1 << 0),
  STATE_RCPD          = (1 << 1),
  STATE_STABLE        = (1 << 2),
  STATE_STABLE_A      = (1 << 3),
  STATE_NYSTROM       = (1 << 4),
  STATE_LOCAL         = (1 << 5),
  STATE_END           = (1 << 6)
};

int dposv_ (char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
int dgesdd_(char *jobz, int *M,int *N,double* A, int* lda, double *S, double *U, int*ldu,
            double* Vt, int *ldvt, double *work, int *lwork, int *iwork, int *info);

static inline void addelem (double *a, int N, double e){ int i; for(i=N-2;i>=0;i--){a[i+1]=a[i];} a[0]=e;}

void rigid(double *H, double *u, const double *G, const double r, dldsize size, const int mode){
  int i,d,k,m,D,K,M,MD; double c[3]; D=size.D;K=size.K;M=size.M; MD=M*D;
  for(m=0;m<M;m++)for(k=0;k<K;k++){
    for(d=0;d<D;d++){c[d]=0;for(i=0;i<D;i++) c[d]+=G[mode?(d+D*i):(i+D*d)]*H[i+D*m+MD*k];}
    for(d=0;d<D;d++){H[d+D*m+k*M*D]=(mode?r:1/r)*c[d];} if(k) continue;
    for(d=0;d<D;d++){c[d]=0;for(i=0;i<D;i++) c[d]+=G[mode?(d+D*i):(i+D*d)]*u[i+D*m];}
    for(d=0;d<D;d++) u[d+D*m]  =(mode?r:1/r)*c[d];
  }
}


static inline void print_search_stats(int loop,double Np,double sigma,double diff,double dslope,int nsgm,dldprms prms,nystmode mode,int state){
  char *str; int direct=!(prms.nyst>0||prms.err>0),ifgt=prms.err>0?1:0;
  str=direct?"direct":(ifgt?"ifgt":(mode&NYST_MODE_SAMPLING?"nystrom":"neighbor (in nystrom)"));

  if(prms.opt&OPT_HISTORY){
    if(loop) printf("   loop=%.3d  mode=%s  sigma=%lf  diff=%lf\n",loop+1,str,sigma,diff);
    else     printf("   loop=%.3d  mode=%s  sigma=%lf          \n",loop+1,str,sigma);
  }
  else{
    if(loop) printf("\033[7A");
    printf("  EM: %d loops\n",    loop+1);
    printf("    Np     = %lf\n",  Np);
    printf("    sigma  = %lf\n",  sigma);
    printf("    conv   = %lf\n",  diff);
    if(loop<nsgm) printf("    dslope = NA\n");
    else          printf("    dslope = %lf\n",dslope);
    printf("    mode   = %s\n", str);
    switch(state){
      case(STATE_START   ): printf("    state  = start   \n"); break;
      case(STATE_RCPD    ): printf("    state  = rcpd    \n"); break;
      case(STATE_NYSTROM ): printf("    state  = nystrom \n"); break;
      case(STATE_LOCAL   ): printf("    state  = local   \n"); break;
      case(STATE_STABLE  ): printf("    state  = stable  \n"); break;
      case(STATE_STABLE_A): printf("    state  = stable+ \n"); break;
      case(STATE_END     ): printf("    state  = end     \n"); break;
    }
  }

  return;
}

double volume(const double *x, int D, int N){
  int d,n; double max=0,min=0,V=1.0;
  for(d=0;d<D;d++){
    for(n=0;n<N;n++) max=n==0?x[d]:fmax(max,x[d+D*n]);
    for(n=0;n<N;n++) min=n==0?x[d]:fmin(min,x[d+D*n]);
    V*=(max-min)*(N+1)/(double)N;
  }
  return V;
}

void state_checker(int state, int flags){
  int on=0; int off=1;
  assert(!((flags&FLAG_IFGT) &&(flags&FLAG_NYSTROM)));
  assert(!((flags&FLAG_IFGT) &&(flags&FLAG_GEMAN  )));
  assert(!((flags&FLAG_GEMAN)&&(flags&FLAG_NYSTROM)));

  switch(state){
    case(STATE_START   ): on=1; off=0; break;
    case(STATE_RCPD    ): on=flags& FLAG_SIMILAR;              off=flags&(FLAG_DEFORM |FLAG_LOCAL|FLAG_ADAPT); break;
    case(STATE_LOCAL   ): on=flags&(FLAG_DEFORM|FLAG_LOCAL);   off=flags&(FLAG_NYSTROM|FLAG_IFGT);             break;
    case(STATE_STABLE  ): on=flags& FLAG_DEFORM;               off=flags&(FLAG_NYSTROM|FLAG_LOCAL|FLAG_ADAPT); break;
    case(STATE_STABLE_A): on=flags&(FLAG_DEFORM|FLAG_ADAPT);   off=flags&(FLAG_NYSTROM|FLAG_LOCAL);            break;
    case(STATE_NYSTROM ): on=flags&(FLAG_DEFORM|FLAG_NYSTROM); off=flags&(FLAG_LOCAL|FLAG_IFGT|FLAG_ADAPT);    break;
    case(STATE_END     ): on=flags& FLAG_DEFORM ;              off=flags& FLAG_NYSTROM;                        break;
    default: printf("Error in state_checker (state=%d). Abort.\n",state); exit(EXIT_FAILURE); break;
  } assert(on); assert(!off);

}

void state_handler(int *state, int *flags, int *flat, int *conv, int *local, dldopts opt){
  int IFOPT_NYST  =  (opt&OPT_NYSTROM  ) ? FLAG_NYSTROM : 0;
  int IFOPT_IFGT  =  (opt&OPT_IFGT     ) ? FLAG_IFGT    : 0;
  int IFOPT_GEMAN =  (opt&OPT_GEMAN    ) ? FLAG_GEMAN   : 0;
  int IFOPT_BARE  = !(opt&OPT_BAREBONES) ? FLAG_SIMILAR : 0;
  int IFOPT_ADAPT =  (opt&OPT_ADAPT    ) ? FLAG_ADAPT   : 0;

  switch(*state){
    case(STATE_START): *conv=0;*local=0;*flat=0;
      if      (opt&OPT_PREALIGN) {*state=STATE_RCPD;    *flags=FLAG_SIMILAR|IFOPT_NYST|IFOPT_IFGT;}
      else if (opt&OPT_NYSTROM ) {*state=STATE_NYSTROM; *flags=FLAG_DEFORM |FLAG_NYSTROM|IFOPT_BARE;}
      else                       {*state=STATE_STABLE;  *flags=FLAG_DEFORM |IFOPT_IFGT|IFOPT_BARE|IFOPT_GEMAN;}
      break;
    case(STATE_RCPD):    if(!(*flat)||(*conv)){break;}  *flat=0;*conv=0;*local=0;
      if      (opt&OPT_NYSTROM ) {*state=STATE_NYSTROM; *flags=FLAG_DEFORM |FLAG_NYSTROM|IFOPT_BARE;}
      else                       {*state=STATE_STABLE;  *flags=FLAG_DEFORM |IFOPT_IFGT|IFOPT_BARE|IFOPT_GEMAN;}
      break;
    case(STATE_NYSTROM): if(!(*flat||*local)){break;}   *flat=0;*conv=0;*local=1;/*CAUTION*/
      *state=STATE_LOCAL; *flags&=~FLAG_NYSTROM; *flags|=FLAG_LOCAL|IFOPT_ADAPT;
      break;
    case(STATE_LOCAL):   if(!(*conv)){break;}
      *state=STATE_END;
      break;
    case(STATE_STABLE):  if(!(*flat||*local)){break;}   *flat=0;*local=0;
      if      (opt&OPT_ADAPT)    {*state=STATE_STABLE_A;*flags|=FLAG_ADAPT;}
      else if (*conv)            {*state=STATE_END;}
      break;
    case(STATE_STABLE_A):if(!(*conv)){break;} *state=STATE_END;
      break;
  }
}

void dld(
    double        *   y,      /*  O  | D x M                | transformed point set     */
    double        *   sgm2,   /*  O  |   1                  | residual sum of squares   */
    double        *   Np,     /*  O  |   1                  | effective num of matching */
    double        *   H,      /* I/O | D x M x K            | model: H (variation)      */
    double        *   u,      /* I/O | D x M                | model: u (mean shape)     */
    dldwork       *   w,      /*  W  | various sizes        | working memory            */
    const double  *   x,      /*  I  | D x N                | target point set          */
    dldsize           size,   /*  I  | N, D, K, M           | matrix size               */
    dldprms           prms    /*  I  | omg,gma,cnv,nlp,opts | dld parameters            */
  ){
  int    i,j,k,m,n,d,M,N,D,K,lp,info,one=1,MD; char uplo='U'; char jobz='A'; int *wi,ws;
  double *Q,*z,*S,*pH,*px,*pu,*t; double *Yc,*Xc,*PXc; double *mX,*mY,*A,*B,*R,*a,*L,*U,*Vt,*wd;
  double c1,c2,s; double diff,noise,val,pres1=1e10,pres2=1e20; double val1,val2,vol;
  double omg,gma,cnv,gem; int nlp,verb,save; int state=STATE_START; int flags=0; int mode=0;
  double *v,*p,*q,*P,*P1,*Pt1,*PX,*xp,*up; ifgtwork *iw; nystwork *nw; double h,e=prms.err;
  int ns=prms.nyst,local,flat,conv; double dlt=4.0,lgma=1e-7; int nsgm=10; double dslope,sgms[256]={0};

  state_handler(&state,&flags,&flat,&conv,&local,prms.opt);

  /* rename variables */
  N=size.N; D=size.D; omg=prms.omg; cnv=prms.cnv; verb=prms.opt & OPT_VERBOSE;  gem=prms.gem;
  K=size.K; M=size.M; gma=prms.gma; nlp=prms.nlp; save=prms.opt & OPT_SAVEPATH;
  MD=M*D;

  Q  =w->Q;  pu=w->pu; PXc=w->PXc; Yc=w->Yc; A=w->A; a=w->a; Vt=w->Vt; S =w->S;
  pH =w->pH; px=w->px; Pt1=w->Pt1; Xc=w->Xc; B=w->B; t=w->t; ws=w->ws; P =w->P; nw=w->nw;
  Yc =w->Yc; xp=w->xp; P1 =w->P1;  mY=w->mY; R=w->R; L=w->L; wi=w->wi; q =w->q; v =w->v;
  Xc =w->Xc; up=w->up; PX =w->PX;  mX=w->mX; U=w->U; z=w->z; wd=w->wd; p =w->p; iw=w->iw;

  /* initialization */
  *sgm2=0;
  for(m=0;m<M;m++)for(d=0;d<D;d++) y[d+D*m]=u[d+D*m];
  for(d=0;d<D;d++){ val1=val2=0;
    for(m=0;m<M;m++)*sgm2+=N*SQ(y[d+D*m]);
    for(n=0;n<N;n++)*sgm2+=M*SQ(x[d+D*n]);
    for(m=0;m<M;m++) val1+=y[d+D*m];
    for(n=0;n<N;n++) val2+=x[d+D*n];
    *sgm2-=2*val1*val2;
  } *sgm2/=(size_t)M*N*D; addelem(sgms,nsgm,sqrt(*sgm2));
  vol=volume(x,D,N);

  /* main computation */
  for(lp=0;lp<nlp;lp++){noise=(pow(2.0*M_PI*(*sgm2),0.5*D)*M*omg)/(vol*(1-omg));h=sqrt(2*(*sgm2));*Np=0;
    if(save)for(m=0;m<M;m++)for(d=0;d<D;d++) S[d+D*m+lp*M*D]=y[d+D*m];
    gma=(flags&FLAG_ADAPT)?lgma:prms.gma;

    if((flags&FLAG_LOCAL)||(flags&FLAG_NYSTROM)||(flags&FLAG_IFGT)){
      /*---------------------------------------------------------------o
      |   approximated computation of P1, Pt1, PX                      |
      o---------------------------------------------------------------*/
      mode=(flags&FLAG_LOCAL)?NYST_MODE_NEIGHBOR:NYST_MODE_SAMPLING;
      for(m=0;m<M;m++) p[m]=1.0;
      if(ns) nystrom(q, nw,x,y,p,D,N,M,ns,h,dlt,mode); else ifgt(q, iw,x,y,p,D,N,M,KLIM,PLIM,h,e);
      for(n=0;n<N;n++) q[n]  =1.0/(q[n]+noise);
      for(n=0;n<N;n++) Pt1[n]=1.0-noise*q[n];
      mode|=NYST_FLAG_REUSE;
      if(ns) nystrom(P1,nw,y,x,q,D,M,N,ns,h,dlt,mode); else ifgt(P1,iw,y,x,q,D,M,N,KLIM,PLIM,h,e);
      for(d=0;d<D;d++){
        for(n=0;n<N;n++) v[n]=q[n]*x[d+D*n];
        if(ns) nystrom(PX+M*d,nw,y,x,v,D,M,N,ns,h,dlt,mode); else ifgt(PX+M*d,iw,y,x,v,D,M,N,KLIM,PLIM,h,e);
      }
      for(m=0;m<M;m++) (*Np)+=P1[m];
    }
    else if(flags&FLAG_GEMAN){
      /*---------------------------------------------------------------o
      |   computing P1, Pt1, PX based on Geman-McClure estimation      |
      o---------------------------------------------------------------*/
      for(n=0;n<N;n++) Pt1[n]=0;
      for(m=0;m<M;m++) P1 [m]=0;
      for(m=0;m<M;m++)for(n=0;n<N;n++){P[m+M*n]=0;for(d=0;d<D;d++)P[m+M*n]+=SQ(x[d+D*n]-y[d+D*m]);}
      for(m=0;m<M;m++)for(n=0;n<N;n++) P[m+M*n]=SQ(gem/(gem+P[m+M*n]));
      for(n=0;n<N;n++)for(m=0;m<M;m++) Pt1[n]+=P[m+M*n];
      for(m=0;m<M;m++)for(n=0;n<N;n++) P1 [m]+=P[m+M*n];
      for(m=0;m<M;m++)for(d=0;d<D;d++){PX[m+M*d]=0;for(n=0;n<N;n++) PX[m+M*d]+=P[m+M*n]*x[d+D*n];}
      for(m=0;m<M;m++) (*Np)+=P1[m];
    }
    else {
      /*---------------------------------------------------------------o
      |   computation of P1, Pt1, PX in direct computation             |
      o---------------------------------------------------------------*/
      for(n=0;n<N;n++) Pt1[n]=0;
      for(m=0;m<M;m++) P1 [m]=0;
      for(m=0;m<M;m++)for(n=0;n<N;n++){val=0;for(d=0;d<D;d++){val+=SQ(x[d+D*n]-y[d+D*m]);} P[m+M*n]=exp(-val/(2.0*(*sgm2)));}
      for(n=0;n<N;n++)for(m=0;m<M;m++) Pt1[n]+=P[m+M*n];
      for(m=0;m<M;m++)for(n=0;n<N;n++) P[m+M*n]/=Pt1[n]+noise;
      for(n=0;n<N;n++) Pt1[n]/=Pt1[n]+noise;
      for(m=0;m<M;m++)for(n=0;n<N;n++) P1 [m]+=P[m+M*n];
      for(m=0;m<M;m++)for(d=0;d<D;d++){PX[m+M*d]=0;for(n=0;n<N;n++) PX[m+M*d]+=P[m+M*n]*x[d+D*n];}
      for(m=0;m<M;m++) (*Np)+=P1[m];
    }

    /*---------------------------------------------o
    |   estimation: shape parameters               |
    o---------------------------------------------*/
    if(!(flags&FLAG_DEFORM)) goto skip_shape;
    /* compute Q (s.t. Qz=b, of size K x K) */
    for(i=0;i<K;i++)for(j=0;j<K;j++) Q[i+K*j]=i==j?gma:0;
    for(m=0;m<M;m++)for(i=0;i<K;i++)for(j=0;j<K;j++){val=0;for(d=0;d<D;d++)val+=H[d+D*m+MD*i]*H[d+D*m+MD*j];Q[i+K*j]+=val*P1[m];}
    for(d=0;d<D;d++)for(k=0;k<K;k++){val=0;for(m=0;m<M;m++)val+=P1[m]*H[d+D*m+MD*k];pH[d+D*k]=val/(*Np);}
    for(i=0;i<K;i++)for(j=0;j<K;j++){val=0;for(d=0;d<D;d++)val+=pH[d+D*i]*pH[d+D*j];Q[i+K*j]-=val*(*Np);}
    for(i=0;i<K;i++)for(j=0;j<K;j++) Q[i+K*j]=(Q[i+K*j]+Q[j+K*i])/2;

    /* compute b (s.t. Qz=b; stored in z) */
    for(m=0;m<M;m++)for(d=0;d<D;d++) xp[d+D*m]=PX[m+M*d];
    for(m=0;m<M;m++)for(d=0;d<D;d++) up[d+D*m]=P1[m]*u[d+D*m];
    for(k=0;k<K;k++){z [k]=0;for(m=0;m<M;m++)for(d=0;d<D;d++)z[k]+=H[d+D*m+MD*k]*(xp[d+D*m]-up[d+D*m]);}
    for(d=0;d<D;d++){px[d]=0;for(n=0;n<N;n++)px[d]+=Pt1[n]*x[d+D*n]; px[d]/=(*Np);}
    for(d=0;d<D;d++){pu[d]=0;for(m=0;m<M;m++)pu[d]+=P1 [m]*u[d+D*m]; pu[d]/=(*Np);}
    for(k=0;k<K;k++){val=0;for(d=0;d<D;d++){val+=pH[d+D*k]*(px[d]-pu[d]);} z[k]-=val*(*Np);}

    /* solve Q*z=b, i.e., z=inv(Q)*b  */
    dposv_(&uplo,&K,&one,Q,&K,z,&K,&info); if(info) halt_execution(ERR_LAPACK_DPOSV);

    /* compute Y (current deformation) */
    for(d=0;d<D;d++){t[d]=px[d]-pu[d];for(k=0;k<K;k++)t[d]-=pH[d+D*k]*z[k];}
    for(m=0;m<M;m++)for(d=0;d<D;d++){y[d+D*m]=u[d+D*m]+t[d];for(k=0;k<K;k++) y[d+D*m]+=H[d+D*m+MD*k]*z[k];}
    skip_shape:

    /*---------------------------------------------o
    |   estimation: pose parameters                |
    o---------------------------------------------*/
    if(!(flags&FLAG_SIMILAR)) goto skip_pose;
    /* centerize X and Y */
    for(d=0;d<D;d++){mX[d]=0;for(n=0;n<N;n++) mX[d]+=x[d+D*n]*Pt1[n];mX[d]/=(*Np);}
    for(d=0;d<D;d++){mY[d]=0;for(m=0;m<M;m++) mY[d]+=y[d+D*m]*P1 [m];mY[d]/=(*Np);}
    for(n=0;n<N;n++)for(d=0;d<D;d++) Xc[n+N*d]=x[d+D*n]-mX[d];
    for(m=0;m<M;m++)for(d=0;d<D;d++) Yc[m+M*d]=y[d+D*m]-mY[d];

    /* A=Xc'*P'*Yc */
    for(m=0;m<M;m++)for(d=0;d<D;d++) PXc[m+M*d]=PX[m+M*d]-P1[m]*px[d];
    for(d=0;d<D;d++)for(i=0;i<D;i++){A[d+i*D]=0;for(m=0;m<M;m++)A[d+i*D]+=PXc[m+M*d]*Yc[m+M*i];}
    for(d=0;d<D;d++)for(i=0;i<D;i++) B[d+i*D]=A[d+i*D];

    /* compute svd of A and rotation matrix R */
    dgesdd_(&jobz,&D,&D,B,&D,L,U,&D,Vt,&D,wd,&ws,wi,&info); if(info) halt_execution(ERR_LAPACK_DGESDD);
    val=det(U,D)*det(Vt,D); if(val<0)for(d=0;d<D;d++)U[d+D*(D-1)]*=-1;
    for(i=0;i<D;i++)for(j=0;j<D;j++){R[i+D*j]=0;for(d=0;d<D;d++) R[i+D*j]+=U[i+d*D]*Vt[d+j*D];}

    /* compute scaling s and intercept a */
    c1=c2=0;
    for(d=0;d<D;d++)for(i=0;i<D;i++)c1+=A[d+i*D]*R[d+D*i];
    for(d=0;d<D;d++)for(m=0;m<M;m++)c2+=SQ(Yc[m+M*d])*P1[m];
    s=c1/c2;

    /* compute transformation */
    for(d=0;d<D;d++){val=0;for(i=0;i<D;i++)val+=R[d+D*i]*mY[i];a[d]=mX[d]-s*val;}
    for(m=0;m<M;m++)for(d=0;d<D;d++){val=0;for(i=0;i<D;i++)val+=R[d+D*i]*y[i+D*m];y[d+D*m]=s*val+a[d];}
    rigid(H,u,R,s,size,1);
    skip_pose:

    /*---------------------------------------------o
    |   estimation: residiual (sgm2)               |
    o---------------------------------------------*/
    pres2=pres1;pres1=*sgm2;*sgm2=0;
    for(d=0;d<D;d++){ val=0;
      for(m=0;m<M;m++)*sgm2+=SQ(y[d+D*m])*P1 [m];
      for(n=0;n<N;n++)*sgm2+=SQ(x[d+D*n])*Pt1[n];
      for(m=0;m<M;m++) val +=y[d+D*m]*PX[m+M*d];
      *sgm2-=2*val;
    } *sgm2/=(*Np)*D;

    /*---------------------------------------------o
    |   check convergence                          |
    o---------------------------------------------*/
    /* dslope: residual slope */
    addelem(sgms,nsgm,sqrt(*sgm2));
    dslope=(sgms[nsgm-1]-sgms[0])/(nsgm*sqrt(*sgm2));
    diff=log(pres2)-log(*sgm2); *sgm2=*sgm2>1e-15?*sgm2:1e-15;

    /* gates */
    flat  =(lp>nsgm) && (dslope<0.01);
    local =(sqrt(*sgm2)*dlt)<0.20;
    conv  =(flags&FLAG_GEMAN)?0:(fabs(diff)<cnv);

    state_handler(&state,&flags,&flat,&conv,&local,prms.opt);
    state_checker(state,flags);

    if(verb) print_search_stats(lp,*Np,sqrt(*sgm2),fabs(diff),dslope,nsgm,prms,mode,state);
    if(state==STATE_END) break;

  } w->lp=lp;

}

dldwork* dldwork_alloc(dldsize size, dldprms prms){
  dldwork *w; int fgt=prms.err>0?1:0,nys=prms.nyst?1:0,T=prms.nyst;
  int N,D,K,M,nlp,save,ws,si=sizeof(int),sd=sizeof(double);

  N=size.N; D=size.D; nlp =prms.nlp;
  K=size.K; M=size.M; save=prms.opt&OPT_SAVEPATH;

  w=calloc(1,sizeof(dldwork)); ws=w->ws=10*D*D;
  w->Q  =calloc(K*K,sd); w->A =calloc(D*D,sd); w->z  =calloc(K,sd); w->a =calloc(D, sd);
  w->pH =calloc(D*K,sd); w->B =calloc(D*D,sd); w->px =calloc(D,sd); w->t =calloc(D, sd);
  w->Yc =calloc(M*D,sd); w->R =calloc(D*D,sd); w->pu =calloc(D,sd); w->L =calloc(D, sd);
  w->Xc =calloc(N*D,sd); w->U =calloc(D*D,sd); w->mX =calloc(D,sd); w->wd=calloc(ws,sd);
  w->PXc=calloc(M*D,sd); w->Vt=calloc(D*D,sd); w->mY =calloc(D,sd); w->wi=calloc(10*D,si);

  w->P1 =calloc(M,sd); w->xp=calloc(M*D,sd); w->PX=calloc(M*D,sd);
  w->Pt1=calloc(N,sd); w->up=calloc(M*D,sd);

  w->p =calloc(M,sd); w->nw =(!nys)?NULL:nystwork_alloc(D,M,N,T);
  w->q =calloc(N,sd); w->iw =(!fgt)?NULL:ifgtwork_alloc(M,N,D,KLIM,PLIM);
  w->v =calloc(N,sd); w->P  =(fgt||nys)?NULL:calloc((size_t)M*N,sd);

  w->S  =save?calloc((size_t)nlp*M*D,sd):NULL;

  return w;
}

dldargs* dldargs_alloc(dldsize size, dldprms prms){
  int N,D,K,M; dldargs *a; double sd=sizeof(double);

  N=size.N; D=size.D;
  K=size.K; M=size.M;

  a=calloc(1,sizeof(dldargs));
  a->y =calloc(D*M,sd);
  a->u =calloc(D*M,sd); a->H=calloc(D*M*K,sd);
  a->x =calloc(D*N,sd); a->w=dldwork_alloc(size,prms);

  return a;
}

void dldargs_input(dldargs *a, const double *x, const double *H, const double *u, dldsize size, dldprms prms){
  int n,d,k,m,N,D,K,M,MD;
  N=size.N; D=size.D; a->prms=prms;
  K=size.K; M=size.M; a->size=size; MD=M*D;
  for(n=0;n<N;n++)for(d=0;d<D;d++) a->x[d+D*n]=x[d+D*n];
  for(m=0;m<M;m++)for(d=0;d<D;d++) a->u[d+D*m]=u[d+D*m];
  for(m=0;m<M;m++)for(d=0;d<D;d++)for(k=0;k<K;k++) a->H[d+D*m+MD*k]=H[d+D*m+MD*k];
}

double dld_wrapper(double *q, int D, void *work, const void *args){
  double Np,sgm2,r,G[9]={0},z[3]={1,0,0}; dldargs *a=(dldargs*) work;
  double *H=a->H,*u=a->u; dldsize size=a->size; dldprms prms=a->prms;
  assert(D==a->size.D);

  if     (D==2) r=alignmat2(G,z,q);
  else if(D==3) r=alignmat3(G,z,q);
  else {printf("D in 'dld_wrapper' must be 2 or 3. Abort.\n");exit(1);}

  #define CD const double
  rigid (H,u,G,r,size,1); /* rigid transformation */
  dld   (a->y,&sgm2,&Np,H,u,a->w,(CD*)a->x,size,prms);
  rigid (H,u,G,r,size,0); /* revert transformation */
  #undef CD

  return Np/sgm2;
}

