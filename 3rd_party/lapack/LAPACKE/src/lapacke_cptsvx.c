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
*****************************************************************************
* Contents: Native high-level C interface to LAPACK function cptsvx
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_cptsvx( int matrix_layout, char fact, lapack_int n,
                           lapack_int nrhs, const float* d,
                           const lapack_complex_float* e, float* df,
                           lapack_complex_float* ef,
                           const lapack_complex_float* b, lapack_int ldb,
                           lapack_complex_float* x, lapack_int ldx,
                           float* rcond, float* ferr, float* berr )
{
    lapack_int info = 0;
    float* rwork = NULL;
    lapack_complex_float* work = NULL;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        LAPACKE_xerbla( "LAPACKE_cptsvx", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        if( LAPACKE_cge_nancheck( matrix_layout, n, nrhs, b, ldb ) ) {
            return -9;
        }
        if( LAPACKE_s_nancheck( n, d, 1 ) ) {
            return -5;
        }
        if( LAPACKE_lsame( fact, 'f' ) ) {
            if( LAPACKE_s_nancheck( n, df, 1 ) ) {
                return -7;
            }
        }
        if( LAPACKE_c_nancheck( n-1, e, 1 ) ) {
            return -6;
        }
        if( LAPACKE_lsame( fact, 'f' ) ) {
            if( LAPACKE_c_nancheck( n-1, ef, 1 ) ) {
                return -8;
            }
        }
    }
#endif
    /* Allocate memory for working array(s) */
    rwork = (float*)LAPACKE_malloc( sizeof(float) * MAX(1,n) );
    if( rwork == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    work = (lapack_complex_float*)
        LAPACKE_malloc( sizeof(lapack_complex_float) * MAX(1,n) );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_1;
    }
    /* Call middle-level interface */
    info = LAPACKE_cptsvx_work( matrix_layout, fact, n, nrhs, d, e, df, ef, b,
                                ldb, x, ldx, rcond, ferr, berr, work, rwork );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_1:
    LAPACKE_free( rwork );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        LAPACKE_xerbla( "LAPACKE_cptsvx", info );
    }
    return info;
}
