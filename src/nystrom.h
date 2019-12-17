typedef struct nystwork{ 
  double *Z/*LD*/; int *U/*L*/; double *v/*P*/; double *K/*PP*/; sortbox *sb/*L*/;
  int *tree; int *ibuf; double *dbuf; int *a; int *t; int **S;
} nystwork;

nystwork *nystwork_alloc(int D, int M, int N, int P);

typedef enum nystmode{
  NYST_MODE_NEIGHBOR  =1,
  NYST_MODE_SAMPLING  =2,
  NYST_FLAG_REUSE     =4 
} nystmode;

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
);
