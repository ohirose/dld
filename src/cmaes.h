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

#define CMA_MAXIMIZE 1
#define CMA_MINIMIZE 0

typedef struct cma_prms {double *w;int mu;double mue;double csgm;double dsgm;double cc;double c1;double cmu;} cma_prms;
typedef struct cma_conf {int epoc; int npop; double tol; int verb; unsigned long seed;} cma_conf;
typedef struct cma_size {int N; int D;} cma_size;

typedef struct cma_stat{
  double *f; sortbox *sb; double *x; double *C;  double *m; double sgm; double *z; double *y; double *wc; 
  double *q; double *yw; double *psgm; double *pc; double *Cr; double *A; double *V; double *L; double *b;
} cma_stat;

cma_prms *cma_prms_alloc (cma_size  sz);
cma_stat *cma_stat_alloc (cma_size  sz);
void      cma_size_init  (cma_size *sz,int D);
void      cma_prms_init  (cma_prms *p, cma_size sz);
void      cma_stat_init  (cma_stat *s, cma_size sz, const double *m, const double sgm);

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
);
