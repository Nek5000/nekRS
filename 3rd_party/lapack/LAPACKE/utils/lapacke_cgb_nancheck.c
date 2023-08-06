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

/* Check a matrix for NaN entries. */

lapack_logical LAPACKE_cgb_nancheck( int matrix_layout, lapack_int m,
                                      lapack_int n, lapack_int kl,
                                      lapack_int ku,
                                      const lapack_complex_float *ab,
                                      lapack_int ldab )
{
    lapack_int i, j;

    if( ab == NULL ) return (lapack_logical) 0;

    if( matrix_layout == LAPACK_COL_MAJOR ) {
        for( j = 0; j < n; j++ ) {
            for( i = MAX( ku-j, 0 ); i < MIN( m+ku-j, kl+ku+1 );
                 i++ ) {
                if( LAPACK_CISNAN( ab[i+(size_t)j*ldab] ) )
                     return (lapack_logical) 1;
            }
        }
    } else if ( matrix_layout == LAPACK_ROW_MAJOR ) {
        for( j = 0; j < n; j++ ) {
            for( i = MAX( ku-j, 0 ); i < MIN( m+ku-j, kl+ku+1 ); i++ ) {
                if( LAPACK_CISNAN( ab[(size_t)i*ldab+j] ) )
                     return (lapack_logical) 1;
            }
        }
    }
    return (lapack_logical) 0;
}
