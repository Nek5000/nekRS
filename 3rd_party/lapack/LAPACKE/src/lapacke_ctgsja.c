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
* Contents: Native high-level C interface to LAPACK function ctgsja
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_ctgsja( int matrix_layout, char jobu, char jobv, char jobq,
                           lapack_int m, lapack_int p, lapack_int n,
                           lapack_int k, lapack_int l, lapack_complex_float* a,
                           lapack_int lda, lapack_complex_float* b,
                           lapack_int ldb, float tola, float tolb, float* alpha,
                           float* beta, lapack_complex_float* u, lapack_int ldu,
                           lapack_complex_float* v, lapack_int ldv,
                           lapack_complex_float* q, lapack_int ldq,
                           lapack_int* ncycle )
{
    lapack_int info = 0;
    lapack_complex_float* work = NULL;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        LAPACKE_xerbla( "LAPACKE_ctgsja", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        if( LAPACKE_cge_nancheck( matrix_layout, m, n, a, lda ) ) {
            return -10;
        }
        if( LAPACKE_cge_nancheck( matrix_layout, p, n, b, ldb ) ) {
            return -12;
        }
        if( LAPACKE_lsame( jobq, 'i' ) || LAPACKE_lsame( jobq, 'q' ) ) {
            if( LAPACKE_cge_nancheck( matrix_layout, n, n, q, ldq ) ) {
                return -22;
            }
        }
        if( LAPACKE_s_nancheck( 1, &tola, 1 ) ) {
            return -14;
        }
        if( LAPACKE_s_nancheck( 1, &tolb, 1 ) ) {
            return -15;
        }
        if( LAPACKE_lsame( jobu, 'i' ) || LAPACKE_lsame( jobu, 'u' ) ) {
            if( LAPACKE_cge_nancheck( matrix_layout, m, m, u, ldu ) ) {
                return -18;
            }
        }
        if( LAPACKE_lsame( jobv, 'i' ) || LAPACKE_lsame( jobv, 'v' ) ) {
            if( LAPACKE_cge_nancheck( matrix_layout, p, p, v, ldv ) ) {
                return -20;
            }
        }
    }
#endif
    /* Allocate memory for working array(s) */
    work = (lapack_complex_float*)
        LAPACKE_malloc( sizeof(lapack_complex_float) * MAX(1,2*n) );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    /* Call middle-level interface */
    info = LAPACKE_ctgsja_work( matrix_layout, jobu, jobv, jobq, m, p, n, k, l,
                                a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu,
                                v, ldv, q, ldq, work, ncycle );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        LAPACKE_xerbla( "LAPACKE_ctgsja", info );
    }
    return info;
}
