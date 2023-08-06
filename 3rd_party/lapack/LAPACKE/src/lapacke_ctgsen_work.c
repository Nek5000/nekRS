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
* Contents: Native middle-level C interface to LAPACK function ctgsen
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_ctgsen_work( int matrix_layout, lapack_int ijob,
                                lapack_logical wantq, lapack_logical wantz,
                                const lapack_logical* select, lapack_int n,
                                lapack_complex_float* a, lapack_int lda,
                                lapack_complex_float* b, lapack_int ldb,
                                lapack_complex_float* alpha,
                                lapack_complex_float* beta,
                                lapack_complex_float* q, lapack_int ldq,
                                lapack_complex_float* z, lapack_int ldz,
                                lapack_int* m, float* pl, float* pr, float* dif,
                                lapack_complex_float* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_ctgsen( &ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb,
                       alpha, beta, q, &ldq, z, &ldz, m, pl, pr, dif, work,
                       &lwork, iwork, &liwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int lda_t = MAX(1,n);
        lapack_int ldb_t = MAX(1,n);
        lapack_int ldq_t = MAX(1,n);
        lapack_int ldz_t = MAX(1,n);
        lapack_complex_float* a_t = NULL;
        lapack_complex_float* b_t = NULL;
        lapack_complex_float* q_t = NULL;
        lapack_complex_float* z_t = NULL;
        /* Check leading dimension(s) */
        if( lda < n ) {
            info = -8;
            LAPACKE_xerbla( "LAPACKE_ctgsen_work", info );
            return info;
        }
        if( ldb < n ) {
            info = -10;
            LAPACKE_xerbla( "LAPACKE_ctgsen_work", info );
            return info;
        }
        if( ldq < n ) {
            info = -14;
            LAPACKE_xerbla( "LAPACKE_ctgsen_work", info );
            return info;
        }
        if( ldz < n ) {
            info = -16;
            LAPACKE_xerbla( "LAPACKE_ctgsen_work", info );
            return info;
        }
        /* Query optimal working array(s) size if requested */
        if( liwork == -1 || lwork == -1 ) {
            LAPACK_ctgsen( &ijob, &wantq, &wantz, select, &n, a, &lda_t, b,
                           &ldb_t, alpha, beta, q, &ldq_t, z, &ldz_t, m, pl, pr,
                           dif, work, &lwork, iwork, &liwork, &info );
            return (info < 0) ? (info - 1) : info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * lda_t * MAX(1,n) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        b_t = (lapack_complex_float*)
            LAPACKE_malloc( sizeof(lapack_complex_float) * ldb_t * MAX(1,n) );
        if( b_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        if( wantq ) {
            q_t = (lapack_complex_float*)
                LAPACKE_malloc( sizeof(lapack_complex_float) *
                                ldq_t * MAX(1,n) );
            if( q_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_2;
            }
        }
        if( wantz ) {
            z_t = (lapack_complex_float*)
                LAPACKE_malloc( sizeof(lapack_complex_float) *
                                ldz_t * MAX(1,n) );
            if( z_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_3;
            }
        }
        /* Transpose input matrices */
        LAPACKE_cge_trans( matrix_layout, n, n, a, lda, a_t, lda_t );
        LAPACKE_cge_trans( matrix_layout, n, n, b, ldb, b_t, ldb_t );
        if( wantq ) {
            LAPACKE_cge_trans( matrix_layout, n, n, q, ldq, q_t, ldq_t );
        }
        if( wantz ) {
            LAPACKE_cge_trans( matrix_layout, n, n, z, ldz, z_t, ldz_t );
        }
        /* Call LAPACK function and adjust info */
        LAPACK_ctgsen( &ijob, &wantq, &wantz, select, &n, a_t, &lda_t, b_t,
                       &ldb_t, alpha, beta, q_t, &ldq_t, z_t, &ldz_t, m, pl, pr,
                       dif, work, &lwork, iwork, &liwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, a_t, lda_t, a, lda );
        LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, b_t, ldb_t, b, ldb );
        if( wantq ) {
            LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, q_t, ldq_t, q, ldq );
        }
        if( wantz ) {
            LAPACKE_cge_trans( LAPACK_COL_MAJOR, n, n, z_t, ldz_t, z, ldz );
        }
        /* Release memory and exit */
        if( wantz ) {
            LAPACKE_free( z_t );
        }
exit_level_3:
        if( wantq ) {
            LAPACKE_free( q_t );
        }
exit_level_2:
        LAPACKE_free( b_t );
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_ctgsen_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_ctgsen_work", info );
    }
    return info;
}
