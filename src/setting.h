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

typedef struct setting {
  double omg; double gma; double cnv; double err; int nyst;/* common */
  double tol; int npop; int epoc; int snlp; int seed;      /* cmaes  */
  int nlp; int bare; int rcpd; int rot; int adap;          /* dld    */
  int cma; int lout; int save; int verb; int ver; int hist;/* flags  */
  char fout[256]; char fprm[256];                          /* files  */
  double gem; /* parameter for comparison */
} setting;

void default_setting   (setting *set);
void verify_setting    (setting *set);
void get_options       (setting *set, int argc, char **argv);
void read_setting      (setting *set, const char *file);
void enable_setting    (dldprms *dp, cma_conf *cf, const setting *set, int argc, char **argv);
void print_setting     (setting set);
