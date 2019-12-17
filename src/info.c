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
#include"version.h"
#include"util.h"
#include"ifgt.h"
#include"nystrom.h"
#include"dld.h"
#include"cmaes.h"
#include"setting.h"

void print_usage(void){
  printf("\n"
         "  USAGE: ./dld  <target: X> <model: H> <model: mu> (+ options)              \n\n"
         "  X:  target point set, matrix of size N x D.                               \n"
         "  H:  Shape model (a) variation,  matrix of size MD x K.                    \n"
         "  mu: Shape model (b) mean shape, vector of length MD.                      \n\n"
         "  OPTIONs: Options must be added AFTER the arguments. If the parameter      \n"
         "  file is set as the argument of '-p', other parameters specified by        \n"
         "  command-line options are ignored.                                         \n\n"
         "  ** with arguments                                                         \n"
         "    -w <omega>, -g <gamma>, -n <nloop>, -c <conv>, -o <output file name>    \n"
         "    -p <parameter file name>.                                               \n\n"  
         "  ** without arguments                                                      \n"
         "    -q quiet, -v version, -s save optmization path, -l switch layout        \n"
         "    pattern of the model file 'H' and 'mu'.                                 \n\n"
         "  EXAMPLE: ./dld hand.txt hat_H.txt hat_mu.txt -w 0.9 -g 1e-4 -n 1000       \n\n");
}

void print_version(void){
  printf("\n");
  printf(" DLD version %s (%s).\n",_VERSION_,_DATE_);
  printf(" Copyright (c) Osamu Hirose                                                 \n\n");
}

void print_info(const char **files, dldsize size){
    const char *fX  =files[0]; const int N=size.N;
    const char *fH  =files[1]; const int D=size.D;
    const char *fmu =files[2]; const int K=size.K;
    const char *fout=files[4]; const int M=size.M;

    printf("\n");
    printf("  Input Data:\n");
    printf("    Target point set (matrix):           [%s]\n",            fX   );
    printf("    Shape model (a) variation  (matrix): [%s]\n",            fH   );
    printf("    Shape model (b) mean shape (vector): [%s]\n",            fmu  ); printf("\n");

    printf("  Input Data Information:\n");
    printf("    Size of target point set:                  [%3d,%2d]\n", N,D  );
    printf("    Size of variation matrix:                  [%3d,%2d]\n", M*D,K);
    printf("    Size of mean shape vector:                 [%3d]\n",     M*D  );
    printf("    Dimension of point-set space:              [%3d]\n",     D    );
    printf("    Number of points in mean shape (computed): [%3d]\n\n",   M    );

    printf("  Output: \n");
    printf("    [%s]\n\n",fout);
}

