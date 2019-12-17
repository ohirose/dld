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

int kmeans(
      int           *  z,    /* OUTPUT: cluster id | size N           */
      double        *  m,    /* OUTPUT: means      | size D x K       */
      int           *  l,    /* OUTPUT: #members   | size K           */
      const double  *  X,    /* INPUT : data       | size D x N       */
      const int        N,    /* INPUT : #samples   | size 1           */
      const int        K,    /* INPUT : #clusters  | size 1           */
      const int        D,    /* INPUT : #variables | size 1           */
      const int        nlp   /* INPUT : #loops     | size 1           */
);
