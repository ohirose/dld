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

void mchoose (int *A, int *M, int *S, int N, int K);
int  nCk     (int n, int k);

typedef struct ifgtwork{ 
  /* the arrays defined in the next line may be reused. */
  int *L; int *A; double *B; double *C; int *G; double *U; double *rx; int *T;
  /* these can be forgotten */  
  int **S; int *t; int *a; double *buf; int *g1; int *g2; int *g3;
} ifgtwork;

/* make K and P large enough */
ifgtwork * ifgtwork_alloc(int M, int N, int D, int K, int P);

void ifgt(
       double         *f,     /*  O  |   M   | resulting values              */
       ifgtwork       *w,     /* I/W |       | working memory                */
       const double   *Y,     /*  I  | D x M | input matrix Y                */
       const double   *X,     /*  I  | D x N | input matrix X                */ 
       const double   *q,     /*  I  |   N   | weight vector q               */
       int             D,     /*  I  | const.| dimension                     */
       int             M,     /*  I  | const.| #points in Y                  */
       int             N,     /*  I  | const.| #points in X                  */ 
       int             K,     /*  *  | const.| upper bound of #clusters      */
       int             P,     /*  *  | const.| upper bound of #trancations   */
       double          h,     /*  I  | const.| band width of Gauss func.     */
       double          e      /*  I  | const.| error bound                   */
);
