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
#include"kmeans.h"

#define SQ(x) ((x)*(x))

int kmeans(
      int           *  z,    /* OUTPUT: cluster id | size N           */
      double        *  m,    /* OUTPUT: means      | size D x K       */
      int           *  l,    /* OUTPUT: #members   | size K           */
      const double  *  X,    /* INPUT : data       | size D x N       */
      const int        N,    /* INPUT : #samples   | size 1           */
      const int        K,    /* INPUT : #clusters  | size 1           */
      const int        D,    /* INPUT : #variables | size 1           */
      const int        nlp   /* INPUT : #loops     | size 1           */
  ){

  int    n,d,k,i;
  int    kmin=-1,step=N/K;
  double v,min=0;

  /* initialization */
  for(k=0;k<K;k++)for(d=0;d<D;d++){m[d+D*k]=X[d+D*k*step];assert(k*step+N*d>=0);}

  /* Main computation */
  for(i=0;i<nlp;i++){
    for(n=0;n<N;n++){
      for(k=0;k<K;k++){v=0;
        for(d=0;d<D;d++)v+=SQ(X[d+D*n]-m[d+D*k]);
        if(k==0||v<min){min=v;kmin=k;}
      } 
      z[n]=kmin; 
    }
    for(k=0;k<K;k++)l[k]=0;
    for(n=0;n<N;n++)if(z[n]>=0)l[z[n]]++;
    //for(k=0;k<K;k++)assert(l[k]);

    /* Update means */
    for(k=0;k<K;k++)for(d=0;d<D;d++) m[d+D*k]=0;
    for(n=0;n<N;n++)for(d=0;d<D;d++)if(z[n]>=0) m[d+D*z[n]]+=X[d+D*n];
    for(k=0;k<K;k++)for(d=0;d<D;d++)if(l[k])    m[d+D*k   ]/=(double)l[k];
  }

  return 0;
}

