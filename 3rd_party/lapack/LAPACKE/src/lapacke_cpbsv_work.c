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
* Contents: Native middle-level C interface to LAPACK function cpbsv
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_cpbsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int kd, lapack_int nrhs,
                               lapack_complex_float* ab, lapack_int ldab,
                               lapack_complex_float* b, lapack_int ldb )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_cpbsv( &uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int ldab_t = MAX(1,kd+1);
        lapack_int ldb_t = MAX(1,n);
        lapack_complex_float* ab_t = NULL;
        lapack_complex_float* b_t = NULL;
        /* Check leading dimension(s) */
        if( ldab < n ) {
            info = -7;
            LAPACKE_xerbla( "LAPACKE_cpbsv_work", info );
            return info;
        }
        if( ldb < nrhs ) {
            info = -9;
            LAPACKE_xerbla( "LAPACKE_cpbsv_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        ab_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * ldab_t * MAX(1,n) );
        if( ab_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        b_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) *
                            ldb_t * MAX(1,nrhs) );
        if( b_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        /* Transpose input matrices */
        LAPACKE_cpb_trans( matrix_layout, uplo, n, kd, ab, ldab, ab_t, ldab_t );
        LAPACKE_cge_trans( matrix_layout, n, nrhs, b, ldb, b_t, ldb_t );
        /* Call LAPACK function and adjust info */
        LAPACK_cpbsv( &uplo, &n, &kd, &nrhs, ab_t, &ldab_t, b_t, &ldb_t,
                      &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_cpb_trans( LAPACK_COL_MAJOR, uplo, n, kd, ab_t, ldab_t, ab,
                           ldab );
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, nrhs, b_t, ldb_t, b, ldb );
        /* Release memory and exit */
        LAPACKE_free( b_t );
exit_level_1:
        LAPACKE_free( ab_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_cpbsv_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_cpbsv_work", info );
    }
    return info;
}
