/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

/* compile with C compiler (not C++) */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "gslib.h"
#include "ogs_findpts.h"

#include "ogstypes.h"

void ogsDevFindptsEval_2(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, void *const in, struct findpts_data_2 *const fd,
  const void *const ogs_fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  ogs_findpts_eval_2( out_base,  out_stride,
                     code_base, code_stride,
                     proc_base, proc_stride,
                       el_base,   el_stride,
                        r_base,    r_stride,
                     npt, in, fd, ogs_fd);
}

void ogsDevFindptsEval_3(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, void *const in, struct findpts_data_3 *const fd,
  const void *const ogs_fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  ogs_findpts_eval_3( out_base,  out_stride,
                     code_base, code_stride,
                     proc_base, proc_stride,
                       el_base,   el_stride,
                        r_base,    r_stride,
                     npt, in, fd, ogs_fd);
}
