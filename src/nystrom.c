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
#include"kdtree.h"
#include"util.h"
#include"nystrom.h"

#define SQ(x) ((x)*(x))
#define MAXTREEDEPTH 64

double genrand_res53(void);
void init_genrand(unsigned long s);

int dpotrs_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
int dpotrf_(char *uplo, int *n, double *A, int *lda, int *info);

inline static double dist2(const double *x, const double *y, int D){
  int d; double val=0;
  for(d=0;d<D;d++) val+=SQ(x[d]-y[d]);
  return val;
}

inline static double gauss(const double *x, const double *y, int D, double h2){return exp(-dist2(x,y,D)/h2);}

inline static void sampling(int *idx, sortbox *sb, int L, int K){
  int l,k; assert(L>=K);
  for(l=0;l<L;l++){sb[l].idx=l;sb[l].val=genrand_res53();}
  qsort(sb,L,sizeof(sortbox),cmp_sortbox_a);
  for(k=0;k<K;k++) idx[k]=sb[k].idx;
}

nystwork *nystwork_alloc(int D, int M, int N, int P){
  int L=M+N,sz=N>M?N:M,si=sizeof(int),sd=sizeof(double),ss=sizeof(sortbox); 
  nystwork *w; assert(L>=P);

  w =calloc(1,sizeof(nystwork));
  w->K=calloc(P*P,sd); w->U=calloc(P,si);
  w->Z=calloc(D*L,sd); w->v=calloc(P,sd); w->sb=calloc(L,ss);

  w->tree=calloc(3*sz+1,si);  w->a=calloc  (sz,si);
  w->ibuf=calloc(6*sz,  si);  w->t=calloc  (sz,si);
  w->dbuf=calloc(2*sz,  sd);  w->S=calloc2i(sz,MAXTREEDEPTH); 

  return w;
}

void nystrom(
       double         *f,     /*  O  |   M   | resulting values              */
       nystwork       *w,     /* I/W |       | working memory                */
       const double   *Y,     /*  I  | D x M | input matrix Y                */
       const double   *X,     /*  I  | D x N | input matrix X                */ 
       const double   *q,     /*  I  |   N   | weight vector b               */
       int             D,     /*  I  | const.| dimension                     */
       int             M,     /*  I  | const.| #points in Y                  */
       int             N,     /*  I  | const.| #points in X                  */ 
       int             P,     /*  I  | const.| #nystrom samples              */
       double          h,     /*  I  | const.| band width                    */
       double          dlt,   /*  I  | const.| neighbor width rate for h     */
       nystmode        mode   /*  I  | const.| sampling vs neighbor (+reuse) */
  ){
  int d,m,n,p,i,j; int one=1,L=M+N; int *U,**S,*t,*a; int info;
  double *Z,*K,*v; sortbox *sb; double rad=dlt*h,reg=1e-9,h2=h*h; char uplo='U';
  
  /* working memory */
  Z=w->Z; U=w->U; v=w->v; t =w->t;
  S=w->S; K=w->K; a=w->a; sb=w->sb;

  if(mode&NYST_MODE_NEIGHBOR) goto neighbor;
  if(mode&NYST_FLAG_REUSE   ) goto reuse;

  init_genrand(clock()); assert(mode&NYST_MODE_SAMPLING);
  for(d=0;d<D;d++)for(m=0;m<M;m++) Z[d+D*( m )]=Y[d+D*m];
  for(d=0;d<D;d++)for(n=0;n<N;n++) Z[d+D*(M+n)]=X[d+D*n];
  sampling(U,sb,L,P);
  for(i=0;i<P;i++)for(j=0;j<P;j++) K[i+P*j]=gauss(Z+D*U[i],Z+D*U[j],D,h2)+(i==j?reg:0);
  dpotrf_(&uplo,&P,K,&P,&info); assert(!info); 

  reuse: assert(mode&NYST_MODE_SAMPLING);
  for(p=0;p<P;p++){v[p]=0; for(n=0;n<N;n++)v[p]+=gauss(Z+D*U[p],X+D*n,D,h2)*q[n];}
  dpotrs_(&uplo,&P,&one,K,&P,v,&P,&info); assert(!info);
  for(m=0;m<M;m++){f[m]=0; for(p=0;p<P;p++)f[m]+=gauss(Z+D*U[p],Y+D*m,D,h2)*v[p];}
  return;

  neighbor:
  kdtree(w->tree,w->ibuf,w->dbuf,X,D,N);
  #ifndef OMP_FOR_CMAES
    #pragma omp parallel for 
  #endif
  for(m=0;m<M;m++){a[m]=t[m]=f[m]=0;
    do{
      eballsearch_next(&a[m],S[m],&t[m],Y+D*m,rad,X,w->tree,D,N);
      if(a[m]>=0) f[m]+=q[a[m]]*gauss(X+D*a[m],Y+D*m,D,h2);
    } while(t[m]);
  }

  return;

}

