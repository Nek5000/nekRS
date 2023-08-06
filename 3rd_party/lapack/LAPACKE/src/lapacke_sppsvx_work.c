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
* Contents: Native middle-level C interface to LAPACK function sppsvx
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_sppsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs, float* ap,
                                float* afp, char* equed, float* s, float* b,
                                lapack_int ldb, float* x, lapack_int ldx,
                                float* rcond, float* ferr, float* berr,
                                float* work, lapack_int* iwork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_sppsvx( &fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x,
                       &ldx, rcond, ferr, berr, work, iwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int ldb_t = MAX(1,n);
        lapack_int ldx_t = MAX(1,n);
        float* b_t = NULL;
        float* x_t = NULL;
        float* ap_t = NULL;
        float* afp_t = NULL;
        /* Check leading dimension(s) */
        if( ldb < nrhs ) {
            info = -11;
            LAPACKE_xerbla( "LAPACKE_sppsvx_work", info );
            return info;
        }
        if( ldx < nrhs ) {
            info = -13;
            LAPACKE_xerbla( "LAPACKE_sppsvx_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        b_t = (float*)LAPACKE_malloc( sizeof(float) * ldb_t * MAX(1,nrhs) );
        if( b_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        x_t = (float*)LAPACKE_malloc( sizeof(float) * ldx_t * MAX(1,nrhs) );
        if( x_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        ap_t = (float*)
            LAPACKE_malloc( sizeof(float) * ( MAX(1,n) * MAX(2,n+1) ) / 2 );
        if( ap_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_2;
        }
        afp_t = (float*)
            LAPACKE_malloc( sizeof(float) * ( MAX(1,n) * MAX(2,n+1) ) / 2 );
        if( afp_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_3;
        }
        /* Transpose input matrices */
        LAPACKE_sge_trans( matrix_layout, n, nrhs, b, ldb, b_t, ldb_t );
        LAPACKE_spp_trans( matrix_layout, uplo, n, ap, ap_t );
        if( LAPACKE_lsame( fact, 'f' ) ) {
            LAPACKE_spp_trans( matrix_layout, uplo, n, afp, afp_t );
        }
        /* Call LAPACK function and adjust info */
        LAPACK_sppsvx( &fact, &uplo, &n, &nrhs, ap_t, afp_t, equed, s, b_t,
                       &ldb_t, x_t, &ldx_t, rcond, ferr, berr, work, iwork,
                       &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, nrhs, b_t, ldb_t, b, ldb );
        LAPACKE_sge_trans( LAPACK_COL_MAJOR, n, nrhs, x_t, ldx_t, x, ldx );
        if( LAPACKE_lsame( fact, 'e' ) && LAPACKE_lsame( *equed, 'y' ) ) {
            LAPACKE_spp_trans( LAPACK_COL_MAJOR, uplo, n, ap_t, ap );
        }
        if( LAPACKE_lsame( fact, 'e' ) || LAPACKE_lsame( fact, 'n' ) ) {
            LAPACKE_spp_trans( LAPACK_COL_MAJOR, uplo, n, afp_t, afp );
        }
        /* Release memory and exit */
        LAPACKE_free( afp_t );
exit_level_3:
        LAPACKE_free( ap_t );
exit_level_2:
        LAPACKE_free( x_t );
exit_level_1:
        LAPACKE_free( b_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_sppsvx_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_sppsvx_work", info );
    }
    return info;
}
