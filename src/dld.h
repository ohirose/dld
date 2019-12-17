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

typedef enum dldopts{
  OPT_PREALIGN        = (1 << 0),
  OPT_BAREBONES       = (1 << 1),
  OPT_VERBOSE         = (1 << 2),
  OPT_SAVEPATH        = (1 << 3),
  OPT_LAYOUT          = (1 << 4),
  OPT_NYSTROM         = (1 << 5),
  OPT_IFGT            = (1 << 6),
  OPT_ADAPT           = (1 << 7),
  OPT_GEMAN           = (1 << 8),
  OPT_HISTORY         = (1 << 9)
} dldopts;

typedef struct dldprms{double omg; double gma; double cnv; double err; double gem; int nyst; int nlp; int opt;} dldprms;
typedef struct dldsize{int N; int M; int D; int K;} dldsize;
typedef struct dldwork{
  int lp; double *Q; double *z; double *t; double *pH; double *px; double *pu; double *Yc;
  double *Xc; double *PXc; double *S; double *mX; double *mY; double *A; double *B;
  double *R; double *a; double *L; double *U; double *Vt; double *wd; int* wi; int ws;
  double *P; double *P1; double *Pt1; double *PX; double *xp; double *up;
  double *p; double *q; double *v; ifgtwork *iw; nystwork *nw;
} dldwork;

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
);

typedef struct dldargs{
  double *y; double *H; double  *x; dldsize size;
             double *u; dldwork *w; dldprms prms;
} dldargs;

void     dldargs_input (dldargs *a, const double *x, const double *H, const double *u, dldsize size, dldprms prms);
double   dld_wrapper   (double *q, int D, void *work, const void *args);
dldargs* dldargs_alloc (dldsize size, dldprms prms);

double volume(const double *x, int D, int N);
