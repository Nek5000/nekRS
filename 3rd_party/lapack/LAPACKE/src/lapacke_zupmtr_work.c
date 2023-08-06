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
* Contents: Native middle-level C interface to LAPACK function zupmtr
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_zupmtr_work( int matrix_layout, char side, char uplo,
                                char trans, lapack_int m, lapack_int n,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, lapack_int ldc,
                                lapack_complex_double* work )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_zupmtr( &side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work,
                       &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int r = LAPACKE_lsame( side, 'l' ) ? m : n;
        lapack_int ldc_t = MAX(1,m);
        lapack_complex_double* c_t = NULL;
        lapack_complex_double* ap_t = NULL;
        /* Check leading dimension(s) */
        if( ldc < n ) {
            info = -10;
            LAPACKE_xerbla( "LAPACKE_zupmtr_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        c_t = (lapack_complex_double*)
            LAPACKE_malloc( sizeof(lapack_complex_double) * ldc_t * MAX(1,n) );
        if( c_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        ap_t = (lapack_complex_double*)
            LAPACKE_malloc( sizeof(lapack_complex_double) *
                            ( MAX(1,r) * MAX(2,r+1) ) / 2 );
        if( ap_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        /* Transpose input matrices */
        LAPACKE_zge_trans( matrix_layout, m, n, c, ldc, c_t, ldc_t );
        LAPACKE_zpp_trans( matrix_layout, uplo, r, ap, ap_t );
        /* Call LAPACK function and adjust info */
        LAPACK_zupmtr( &side, &uplo, &trans, &m, &n, ap_t, tau, c_t, &ldc_t,
                       work, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, m, n, c_t, ldc_t, c, ldc );
        /* Release memory and exit */
        LAPACKE_free( ap_t );
exit_level_1:
        LAPACKE_free( c_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_zupmtr_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_zupmtr_work", info );
    }
    return info;
}
