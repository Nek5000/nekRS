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
* Contents: Native middle-level C interface to LAPACK function zgeesx
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_zgeesx_work( int matrix_layout, char jobvs, char sort,
                                LAPACK_Z_SELECT1 select, char sense,
                                lapack_int n, lapack_complex_double* a,
                                lapack_int lda, lapack_int* sdim,
                                lapack_complex_double* w,
                                lapack_complex_double* vs, lapack_int ldvs,
                                double* rconde, double* rcondv,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_logical* bwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_zgeesx( &jobvs, &sort, select, &sense, &n, a, &lda, sdim, w, vs,
                       &ldvs, rconde, rcondv, work, &lwork, rwork, bwork,
                       &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int lda_t = MAX(1,n);
        lapack_int ldvs_t = MAX(1,n);
        lapack_complex_double* a_t = NULL;
        lapack_complex_double* vs_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -8;
            LAPACKE_xerbla( "LAPACKE_zgeesx_work", info );
            return info;
        }
        if( ldvs < n ) {
            info = -12;
            LAPACKE_xerbla( "LAPACKE_zgeesx_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( lwork == -1 ) {
            LAPACK_zgeesx( &jobvs, &sort, select, &sense, &n, a, &lda_t, sdim,
                           w, vs, &ldvs_t, rconde, rcondv, work, &lwork, rwork,
                           bwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (lapack_complex_double*)
            LAPACKE_malloc( sizeof(lapack_complex_double) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        if( LAPACKE_lsame( jobvs, 'v' ) ) {
            vs_t = (lapack_complex_double*)
                LAPACKE_malloc( sizeof(lapack_complex_double) *
                                ldvs_t * MAX(1,n) );
            if( vs_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_1;
            }
        }
        /* Transpose input matrices */
        LAPACKE_zge_trans( matrix_layout, n, n, a, lda, a_t, lda_t );
        /* Call LAPACK function and adjust info */
        LAPACK_zgeesx( &jobvs, &sort, select, &sense, &n, a_t, &lda_t, sdim, w,
                       vs_t, &ldvs_t, rconde, rcondv, work, &lwork, rwork,
                       bwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, a_t, lda_t, a, lda );
        if( LAPACKE_lsame( jobvs, 'v' ) ) {
            LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, n, vs_t, ldvs_t, vs, ldvs );
        }
        /* Release memory and exit */
        if( LAPACKE_lsame( jobvs, 'v' ) ) {
            LAPACKE_free( vs_t );
        }
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_zgeesx_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_zgeesx_work", info );
    }
    return info;
}
