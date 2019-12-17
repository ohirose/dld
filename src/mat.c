// Copyright (c) 2014--2017 Osamu Hirose
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
#include"util.h"

double dot2(const double v1[2], const double v2[2]){
  return v1[0]*v2[0]+v1[1]*v2[1];
}

double norm2(const double v[2]){
  return sqrt(dot2(v,v));
}

void unit2(double v[2]){
  double n=norm2(v); v[0]/=n;v[1]/=n;v[2]/=n;
}

void cross3(double v[3], const double v1[3], const double v2[3]){
  v[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

double dot3(const double v1[3], const double v2[3]){
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double norm3(const double v[3]){
  return sqrt(dot3(v,v));
}

void unit3(double v[3]){
  double n=norm3(v); v[0]/=n;v[1]/=n;v[2]/=n;
}

double dist3(const double v1[3], const double v2[3]){
  double val=0;
  val+=SQ(v1[0]-v2[0]);
  val+=SQ(v1[1]-v2[1]);
  val+=SQ(v1[2]-v2[2]);
  return sqrt(val);
}


double alignmat2(double *G, const double v1[2], const double v2[2]){
  double l1,l2,c,s,u[2],R[2][2],W[2][2],U[2][2]; int i,j;

  for(i=0;i<2;i++)for(j=0;j<2;j++) G[i+2*j]=U[i][j]=R[i][j]=W[i][j]=0;
  /* u: pi/2 rotation of v1 */
  u[0]=-v1[1]; u[1]=v1[0];

  /* rotation matrix */
  c=dot2(v1,v2)/(norm2(v1)*norm2(v2));
  s=sqrt(1-c*c)*(dot2(v2,u)<0?-1.0:1.0);
  R[0][0]=c; R[0][1]=-s; 
  R[1][0]=s; R[1][1]= c;

  /* alignment matrix */
  l1=norm2(v1); l2=norm2(v2);
  for(i=0;i<2;i++)for(j=0;j<2;j++) G[i+2*j]=R[i][j];

  return l2/l1;
}

double alignmat3(double *G, const double v1[3], const double v2[3]){
   double e=1e-20,c,s,w1[3],w2[3],U[3][3],R[3][3],W[3][3]; int i,j,k;

   /* case: either of two vectors is zero */
   if(norm3(v1)<e||norm3(v2)<e){
     printf("Can't align v1 and v2 due to a zero vector in alignmat3. Abort.\n");
     exit(EXIT_FAILURE);
   }

   /* initialization */
   for(i=0;i<3;i++)for(j=0;j<3;j++) G[i+3*j]=U[i][j]=R[i][j]=W[i][j]=0;
   for(i=0;i<3;i++){w1[i]=v1[i];} unit3(w1);
   for(i=0;i<3;i++){w2[i]=v2[i];} unit3(w2);

   /* case: two vectors are parallel */
   if(dist3(w1,w2)<e){G[0]=G[4]=G[8]=1.0; goto final;}

   /* case: otherwise */
   /* basis transformation */
   for(i=0;i<3;i++){U[0][i]=v1[i];} unit3(U[0]);
   cross3(U[2],v1,v2  );            unit3(U[2]);
   cross3(U[1],U[2],v1);            unit3(U[1]);

   /* rotation matrix */
   c=dot3(v1,v2)/(norm3(v1)*norm3(v2));
   s=sqrt(1-c*c)*(dot3(v2,U[1])<0?-1.0:1.0);
   R[0][0]=c; R[0][1]=-s;
   R[1][0]=s; R[1][1]= c; R[2][2]=1;

   /* alignment matrix */
   for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++) W[i][j] +=R[i][k]*U[k][j];
   for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++) G[i+3*j]+=U[k][i]*W[k][j];

   final:
   return norm3(v2)/norm3(v1);
}

double det(const double *A, const int D){
  assert(D==2||D==3);
  return D==2? A[0]*A[3]-A[1]*A[2]:
               A[0+D*0]*(A[1+D*1]*A[2+D*2]-A[1+D*2]*A[2+D*1])
              -A[0+D*1]*(A[1+D*0]*A[2+D*2]-A[1+D*2]*A[2+D*0])
              +A[0+D*2]*(A[1+D*0]*A[2+D*1]-A[1+D*1]*A[2+D*0]);
}

