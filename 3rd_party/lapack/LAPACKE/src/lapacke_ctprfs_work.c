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
* Contents: Native middle-level C interface to LAPACK function ctprfs
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_ctprfs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int nrhs,
                                const lapack_complex_float* ap,
                                const lapack_complex_float* b, lapack_int ldb,
                                const lapack_complex_float* x, lapack_int ldx,
                                float* ferr, float* berr,
                                lapack_complex_float* work, float* rwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_ctprfs( &uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx,
                       ferr, berr, work, rwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int ldb_t = MAX(1,n);
        lapack_int ldx_t = MAX(1,n);
        lapack_complex_float* b_t = NULL;
        lapack_complex_float* x_t = NULL;
        lapack_complex_float* ap_t = NULL;
        /* Check leading dimension(s) */
        if( ldb < nrhs ) {
            info = -9;
            LAPACKE_xerbla( "LAPACKE_ctprfs_work", info );
            return info;
        }
        if( ldx < nrhs ) {
            info = -11;
            LAPACKE_xerbla( "LAPACKE_ctprfs_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        b_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) *
                            ldb_t * MAX(1,nrhs) );
        if( b_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        x_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) *
                            ldx_t * MAX(1,nrhs) );
        if( x_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        ap_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) *
                            ( MAX(1,n) * MAX(2,n+1) ) / 2 );
        if( ap_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_2;
        }
        /* Transpose input matrices */
        LAPACKE_cge_trans( matrix_layout, n, nrhs, b, ldb, b_t, ldb_t );
        LAPACKE_cge_trans( matrix_layout, n, nrhs, x, ldx, x_t, ldx_t );
        LAPACKE_ctp_trans( matrix_layout, uplo, diag, n, ap, ap_t );
        /* Call LAPACK function and adjust info */
        LAPACK_ctprfs( &uplo, &trans, &diag, &n, &nrhs, ap_t, b_t, &ldb_t, x_t,
                       &ldx_t, ferr, berr, work, rwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Release memory and exit */
        LAPACKE_free( ap_t );
exit_level_2:
        LAPACKE_free( x_t );
exit_level_1:
        LAPACKE_free( b_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_ctprfs_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_ctprfs_work", info );
    }
    return info;
}
