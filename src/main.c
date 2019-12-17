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
#include<ctype.h>
#include<math.h>
#include<unistd.h>
#include"util.h"
#include"ifgt.h"
#include"nystrom.h"
#include"dld.h"
#include"cmaes.h"
#include"setting.h"
#include"info.h"
#include"error.h"

#define CD const double
#define MAXLENGTH 10000
#define MAXLEN    1024

int save_opt_path(const char *fpath, const double *x, const dldargs *a, dldsize sz){
  int N=sz.N,M=sz.M,D=sz.D,lp=a->w->lp; double *S=a->w->S; int si=sizeof(int),sd=sizeof(double);
  FILE *fp=fopen(fpath,"wb"); int verb=a->prms.opt&OPT_VERBOSE;
  if(!fp){printf("Can't open: %s\n",fpath);exit(EXIT_FAILURE);}
  fwrite(&N, si,   1,  fp);
  fwrite(&D, si,   1,  fp);
  fwrite(&M, si,   1,  fp);
  fwrite(&lp,si,   1,  fp);
  fwrite(S,  sd,lp*D*M,fp);
  fwrite(x,  sd,   D*N,fp);
  fclose(fp);
  if(verb) printf("\n  ** Search path during optimization was saved to: [%s]\n\n",fpath);
  return 0;
}

int main(int argc, char **argv){
  char mode; int l,m,n,k,d,N,M,K,D,L,num[6]; double *x,*X,*u,*U,*H,*H0,*T;
  char *na="NA",*files[5],ftxt[MAXLEN]="T.txt",fbin[MAXLEN]="T.bin"; setting set;
  dldprms dp[3]; dldsize ds; int cma=0,save=0,lout=0,verb=1; int sod=sizeof(double);
  cma_stat *cst; cma_prms *cpr; cma_conf cf; cma_size csz; dldargs *a,*a0,*a2,**a1;
  double isgm,*im,*q,*q1,*q0; double V,scl; char *fpath=".optpath.bin";

  if(argc==2){print_version();exit(EXIT_SUCCESS);}
  if(argc <4){print_usage();  exit(EXIT_SUCCESS);}

  /* read input files */
  X =read2d_m(num+0,num+1,&mode,argv[1],na);
  H0=read2d_m(num+2,num+3,&mode,argv[2],na);
  U =read2d_m(num+4,num+5,&mode,argv[3],na);
  if(num[1]>3)         halt_execution(ERR_UNSUPPORTED_DIM);
  if(num[2]!=num[4])   halt_execution(ERR_INCOSISTENT_NROW);
  if(num[5]!=1)        halt_execution(ERR_INCOSISTENT_FORMAT);
  if(num[2]%num[1]!=0) halt_execution(ERR_INCOSISTENT_DIM);

  /* set parameters; command-line options are overwritten if parameter file is available. */
  default_setting(&set); get_options (&set,argc,argv);
  if( strlen(set.fprm))  read_setting(&set,set.fprm);
  if(!strlen(set.fout))  strcpy(set.fout,(mode=='t')?ftxt:fbin);
  verify_setting (&set);
  enable_setting (dp,&cf,&set,argc,argv);

  if(set.ver) print_version();

  cma =set.cma;  lout=set.lout;
  verb=set.verb; save=set.save;

  files[0]=argv[1];  N=ds.N=num[0];
  files[1]=argv[2];  D=ds.D=num[1];
  files[2]=argv[3];  K=ds.K=num[3];
  files[3]=set.fprm; M=ds.M=num[2]/D;
  files[4]=set.fout;

  if(verb) print_info    ((const char**)files,ds);
  if(verb) print_setting (set);

  u=calloc(D*M,  sod); q0=calloc(D,sod); q0[0]=1.0;
  x=calloc(D*N,  sod); q1=calloc(D,sod); q=cma?q1:q0;
  H=calloc(D*M*K,sod); im=calloc(D,sod); T=calloc(M*D,sod);

  for(m=0;m<M;m++)for(k=0;k<K;k++)for(d=0;d<D;d++) H[d+D*m+M*D*k]=H0[(lout?(d+D*m):(m+M*d))+M*D*k];
  for(m=0;m<M;m++)for(d=0;d<D;d++) u[d+D*m]=U[lout?(d+D*m):(m+M*d)];
  for(n=0;n<N;n++)for(d=0;d<D;d++) x[d+D*n]=X[n+N*d];

  V=volume(u,D,M); scl=pow(V,1.0/(double)D);
  for(m=0;m<M;m++)for(k=0;k<K;k++)for(d=0;d<D;d++) H[d+D*m+M*D*k]/=scl;
  for(m=0;m<M;m++)for(d=0;d<D;d++) u[d+D*m]/=scl;
  for(n=0;n<N;n++)for(d=0;d<D;d++) x[d+D*n]/=scl;

  /* dld memory allocation */
  a0=dldargs_alloc(ds,dp[0]); dldargs_input(a0,(CD*)x,(CD*)H,(CD*)u,ds,dp[0]);
  a2=dldargs_alloc(ds,dp[2]); dldargs_input(a2,(CD*)x,(CD*)H,(CD*)u,ds,dp[2]);
  a =cma?a2:a0;

  if(!cma) goto final;

  /* cma initialization & allocation */
  cma_size_init(&csz,D);   isgm=0.1; L=csz.N=set.npop?cf.npop:csz.N; //im[0]=1.0;
  cpr=cma_prms_alloc(csz); cma_prms_init(cpr,csz);
  cst=cma_stat_alloc(csz); cma_stat_init(cst,csz,im,isgm);
  a1 =calloc(L,sizeof(dldargs**));
  for(l=0;l<L;l++) a1[l]=dldargs_alloc(ds,dp[1]);
  for(l=0;l<L;l++) dldargs_input(a1[l],(CD*)x,(CD*)H,(CD*)u,ds,dp[1]);

  /* cma main computation */
  cma_optimize(q,cst,(void**)a1,NULL,dld_wrapper,cpr,cf,csz,CMA_MAXIMIZE);

  /* dld main computation */
  final: dld_wrapper (q,D,(void*)a,NULL);

  /* file output */
  if(save) save_opt_path(fpath,x,a,ds);
  for(m=0;m<M;m++)for(d=0;d<D;d++) T[m+M*d]=scl*a->y[d+D*m];
  write2d_m(set.fout,(CD*)T,M,D,na);

  return 0;
}

