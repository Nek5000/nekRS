/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************
* Contents: Native C interface to LAPACK utility function
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

/* Converts input general band matrix from row-major(C) to
 * column-major(Fortran) layout or vice versa.
 */

void LAPACKE_zgb_trans( int matrix_layout, lapack_int m, lapack_int n,
                        lapack_int kl, lapack_int ku,
                        const lapack_complex_double *in, lapack_int ldin,
                        lapack_complex_double *out, lapack_int ldout )
{
    lapack_int i, j;

    if( in == NULL || out == NULL ) return;

    if( matrix_layout == LAPACK_COL_MAJOR ) {
        for( j = 0; j < MIN( ldout, n ); j++ ) {
            for( i = MAX( ku-j, 0 ); i < MIN3( ldin, m+ku-j, kl+ku+1 );
                 i++ ) {
                out[(size_t)i*ldout+j] = in[i+(size_t)j*ldin];
            }
        }
    } else if ( matrix_layout == LAPACK_ROW_MAJOR ) {
        /* TODO: interchange loops for performance.
         * This is just reference impemeltation
         */
        for( j = 0; j < MIN( n, ldin ); j++ ) {
            for( i = MAX( ku-j, 0 ); i < MIN3( ldout, m+ku-j, kl+ku+1 );
                 i++ ) {
                out[i+(size_t)j*ldout] = in[(size_t)i*ldin+j];
            }
        }
    }
}
