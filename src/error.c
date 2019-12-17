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
#include"error.h"

void halt_execution(errtype err){
  char message[6][256]={
    "ERROR: The number of rows for <model: H> and <model: mu> are inconsistent.",
    "ERROR: File format of the mean shape <model: mu> must be a vector of length M x D.",
    "ERROR: The number of rows for <model: H> is not a factor of D.",
    "ERROR: Dimemnsion of the space in which point sets are embedded must be 2 or 3, currently.",
    "ERROR: Eigendecomposition failed. Increasing #samples for Nystrom might help if turned on.",
    "ERROR: SVD failed. This might be caused in a situation that input 3D point sets are embedded into 2D space."
  };

  printf("\n %s Abort.\n\n",message[err]);
  exit(EXIT_FAILURE);
}

