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
#include<unistd.h>
#include<assert.h>
#include"util.h"
#include"ifgt.h"
#include"nystrom.h"
#include"dld.h"
#include"cmaes.h"
#include"setting.h"

#define MAXLENGTH 1024

void default_setting(setting *set){
  set->omg =0.10; set->epoc=300;  set->cma =0; set->lout=0;
  set->gma =1e-3; set->snlp=100;  set->seed=1; set->bare=0;
  set->cnv =1e-4; set->nlp =1000; set->rot =0; set->verb=1;
  set->tol =1e-4; set->npop=0;    set->rcpd=0; set->save=0; 
  set->err =0;    set->ver =0;    set->nyst=0; set->adap=0;
  set->gem =0;    set->hist=0;

  strcpy(set->fprm,"");
  strcpy(set->fout,"");
}

void get_options(setting *set, int argc, char **argv){
  int opt,optc=argc-3; char **optv=argv+3;
  while((opt=getopt(optc,optv,"j:y:f:e:t:m:b:x:w:g:c:n:p:o:arqslvzihH"))!=-1){
    switch(opt){
      case 'w': set->omg =atof(optarg);    break; //  0
      case 'g': set->gma =atof(optarg);    break; //  1
      case 'c': set->cnv =atof(optarg);    break; //  2
      case 't': set->tol =atof(optarg);    break; //  3
      case 'f': set->err =atof(optarg);    break; //  4
      case 'j': set->gem =atof(optarg);    break; //  5
      case 'y': set->nyst=atoi(optarg);    break; //  6
      case 'b': set->npop=atoi(optarg);    break; //  7
      case 'e': set->epoc=atoi(optarg);    break; //  8
      case 'm': set->snlp=atoi(optarg);    break; //  9
      case 'x': set->seed=atoi(optarg);    break; // 10
      case 'n': set->nlp =atoi(optarg);    break; // 11
      case 'r': set->rot =1;               break; // 12
      case 'z': set->bare=1;               break; // 13
      case 'i': set->rcpd=1;               break; // 14
      case 'a': set->cma =1;               break; // 15
      case 'l': set->lout=1;               break; // 16
      case 's': set->save=1;               break; // 17
      case 'q': set->verb=0;               break; // 18
      case 'v': set->ver =1;               break; // 19
      case 'h': set->adap=1;               break; // 20
      case 'H': set->hist=1;               break; // 21
      case 'p': strcpy (set->fprm,optarg); break; // 22
      case 'o': strcpy (set->fout,optarg); break; // 23
      default : exit(EXIT_FAILURE);
    }
  }
}

void read_setting (setting *set, const char *file){
  int j,n,l,L,N=18; FILE *fp; char *s,*b,line[MAXLENGTH]; char *form[]={"%lf","%d","%s"};
  char *names[]={
    "omega","gamma","conv","tol","error","geman","nyst","npop","epoc","snlp","seed","rot",
    "nloop","bare","rcpd","cmaes","layout","save","verb","fout"
  };
 
  fp=fopen(file,"r");if(!fp){printf("File not found: %s\n",file);exit(EXIT_FAILURE);}
  while(fgets(line,MAXLENGTH,fp)){
    if(strchr("#|+\n",line[0])){continue;} L=strlen(line);
    for(l=0;l<L;l++) line[l]=tolower(line[l]);
    for(n=0;n<N;n++){s=names[n];if(!strstr(line,s)) continue; b=line+strlen(s);
      switch(n){
        case  0: j=sscanf(b,form[0],&(set->omg )); break;
        case  1: j=sscanf(b,form[0],&(set->gma )); break;
        case  2: j=sscanf(b,form[0],&(set->cnv )); break;
        case  3: j=sscanf(b,form[0],&(set->tol )); break;
        case  4: j=sscanf(b,form[0],&(set->err )); break;
        case  5: j=sscanf(b,form[0],&(set->gem )); break;
        case  6: j=sscanf(b,form[1],&(set->nyst)); break;
        case  7: j=sscanf(b,form[1],&(set->npop)); break;
        case  8: j=sscanf(b,form[1],&(set->epoc)); break;
        case  9: j=sscanf(b,form[1],&(set->snlp)); break;
        case 10: j=sscanf(b,form[1],&(set->seed)); break;
        case 11: j=sscanf(b,form[1],&(set->rot )); break;
        case 12: j=sscanf(b,form[1],&(set->nlp )); break;
        case 13: j=sscanf(b,form[1],&(set->bare)); break;
        case 14: j=sscanf(b,form[1],&(set->rcpd)); break;
        case 15: j=sscanf(b,form[1],&(set->cma )); break;
        case 16: j=sscanf(b,form[1],&(set->lout)); break;
        case 17: j=sscanf(b,form[1],&(set->save)); break;
        case 18: j=sscanf(b,form[1],&(set->verb)); break;
        case 19: j=sscanf(b,form[1],&(set->adap)); break;
        case 20: j=sscanf(b,form[2],  set->fout ); break;
      } if(j!=1) printf("Can't read the line: %s",line);
    }
  } fclose(fp);
}

void verify_setting(setting *set){
  if(set->omg <=0.0)  {printf("\n ERROR: Parameter 'omega' must be in (0,1). Abort. \n\n"); goto halt;}
  if(set->omg >=1.0)  {printf("\n ERROR: Parameter 'omega' must be in (0,1). Abort. \n\n"); goto halt;}
  if(set->gma <=0.0)  {printf("\n ERROR: Parameter 'gamma' must be positive. Abort. \n\n"); goto halt;}
  if(set->cnv <=0.0)  {printf("\n ERROR: Parameter 'conv' must be positive.  Abort. \n\n"); goto halt;}
  if(set->nlp <=0)    {printf("\n ERROR: Parameter 'nlp' must be positive.   Abort. \n\n"); goto halt;}
  if(set->err < 0.0)  {printf("\n ERROR: Parameter 'error' must be positive or zero. Abort. \n\n"); goto halt;}
  if(set->gem < 0.0)  {printf("\n ERROR: Parameter 'geman' must be positive or zero. Abort. \n\n"); goto halt;}
  if(set->nyst< 0)    {printf("\n ERROR: Parameter 'nyst' must be positive.  Abort. \n\n"); goto halt;}
  if(set->cma){
    if(set->tol <=0.0){printf("\n ERROR: Parameter 'tol' must be positive.   Abort. \n\n"); goto halt;}
    if(set->npop< 0)  {printf("\n ERROR: Parameter 'npop' must be positive.  Abort. \n\n"); goto halt;}
    if(set->epoc<=0)  {printf("\n ERROR: Parameter 'epoc' must be positive.  Abort. \n\n"); goto halt;}
    if(set->snlp<=0)  {printf("\n ERROR: Parameter 'snlp' must be positive.  Abort. \n\n"); goto halt;}
    if(set->seed<=0)  {printf("\n ERROR: Parameter 'seed' must be positive.  Abort. \n\n"); goto halt;}
  }

  if(set->err>0&&set->cma){ set->err=0;
    printf("\n WARNING: Turnning BOTH ifgt and cma-es options on is currently unsupported. Disabled IFGT.\n");
  }
  if(set->nyst>0&&set->err>0){ set->err=0;
    printf("\n WARNING: Can't turnning BOTH ifgt and nystrom options on. Disabled IFGT.\n");
  }
  if(set->nyst>0&&set->gem>0){ set->nyst=0;
    printf("\n WARNING: Can't turnning BOTH Geman-McClure and Nystrom options on. Disabled Nystrom.\n");
  }
  if(set->err >0&&set->gem>0){ set->err =0;
    printf("\n WARNING: Can't turnning BOTH Geman-McClure and IFGT options on. Disabled IFGT.\n");
  }
  return;

  halt: exit(EXIT_FAILURE);
}

void enable_setting(dldprms *dp, cma_conf *cf, const setting *set, int argc, char **argv){
  int ifgt= set->err>0?1:0, nyst=set->nyst?1:0, gem=set->gem?1:0;
  /* default estimation mode (& initialize options) */
  dp[0].opt=dp[1].opt=dp[2].opt=(ifgt?OPT_IFGT:0)|(nyst?OPT_NYSTROM:0)|(gem?OPT_GEMAN:0);

  /* enable dld setting */
  dp[0].omg=dp[1].omg=dp[2].omg=set->omg;
  dp[0].gma=dp[1].gma=dp[2].gma=set->gma;
  dp[0].cnv=dp[1].cnv=dp[2].cnv=set->cnv;
  dp[0].err=dp[1].err=dp[2].err=set->err;
  dp[0].gem=dp[1].gem=dp[2].gem=set->gem;
  dp[0].nyst=dp[1].nyst=dp[2].nyst=set->nyst;
  dp[0].nlp=dp[2].nlp=set->nlp;            
  dp[1].nlp=set->snlp;
  if(!(set->rot)) dp[1].opt|=OPT_BAREBONES;
  if(set->bare)   dp[0].opt|=OPT_BAREBONES;
  if(set->rcpd)   dp[0].opt|=OPT_PREALIGN;
  if(set->save)  {dp[0].opt|=OPT_SAVEPATH;  dp[2].opt|=OPT_SAVEPATH;}
  if(set->verb)  {dp[0].opt|=OPT_VERBOSE;   dp[2].opt|=OPT_VERBOSE;}
  if(set->adap)  {dp[0].opt|=OPT_ADAPT;     dp[2].opt|=OPT_ADAPT;}
  if(set->hist)  {dp[0].opt|=OPT_HISTORY;   dp[2].opt|=OPT_HISTORY;}

  /* enable cmaes configuration */
  cf->npop=set->npop; cf->tol =set->tol;  cf->verb=set->verb;
  cf->epoc=set->epoc; cf->seed=set->seed;
}

void print_setting(setting set){
  printf("  Defomation Mode:\n    ");
  if(set.cma) printf("[Pose: CMAES] + [Shape: DLD]\n\n");
  else {
    if( set.rcpd && !set.bare) printf("[Init: Rigid CPD] + [Shape+Pose: DLD] \n\n");
    if( set.rcpd &&  set.bare) printf("[Init: Rigid CPD] + [Shape: DLD] \n\n");
    if(!set.rcpd && !set.bare) printf("[Shape+Pose: DLD] \n\n");
    if(!set.rcpd &&  set.bare) printf("[Shape: DLD] \n\n");
  }

  if(set.gem>0) printf("  Estimator:\n    [Geman-McClure: parameter = %lf]\n\n",set.gem);
  else          printf("  Estimator:\n    [Maximum Likelihood]\n\n");

  printf("  Computation Mode:\n    "); assert(!(set.nyst>0&&set.err));
  if      (set.nyst)  printf("[Nystrom: #sample = %d]", set.nyst);
  else if (set.err>0) printf("[IFGT: error bound = %E]",set.err);
  else                printf("[Direct]");
  if      (set.adap)  printf(" + [Adaptive gamma control]");
  printf("\n\n");

  printf("  DLD parameters: "); if(strlen(set.fprm))printf("[%s]",set.fprm ); printf("\n");
  printf("    layout     = %s\n",  set.lout?"xyzxyz":"xxyyzz");
  printf("    omega      ="); if(set.omg>1e-4) printf(" %.4lf\n",set.omg); else printf(" %.2E\n",set.omg);
  printf("    gamma      ="); if(set.gma>1e-4) printf(" %.4lf\n",set.gma); else printf(" %.2E\n",set.gma);
  printf("    conv       ="); if(set.cnv>1e-4) printf(" %.4lf\n",set.cnv); else printf(" %.2E\n",set.cnv);
  printf("    nloop      = %d\n",  set.nlp);
  if(!set.cma) goto skip;
  printf("\n");
  printf("  CMAES parameters: \n");
  printf("    epoc       = %d \n", set.epoc);
  printf("    tol        = %lf\n", set.tol);
  printf("    npop       = %d\n",  set.npop);
  printf("    snlp       = %d\n",  set.snlp);
  printf("    seed       = %d\n",  set.seed);
  skip:
  printf("\n");

  if(strlen(set.fprm)){printf("\n");
    printf("  ** Parameters specified by %s were used, that is,\n", set.fprm);
    printf("  those specified by command-line arguments are disabled.\n\n");
  }
}  

