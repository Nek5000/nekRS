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
* Contents: Native middle-level C interface to LAPACK function zhpgvx
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_zhpgvx_work( int matrix_layout, lapack_int itype, char jobz,
                                char range, char uplo, lapack_int n,
                                lapack_complex_double* ap,
                                lapack_complex_double* bp, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_complex_double* work, double* rwork,
                                lapack_int* iwork, lapack_int* ifail )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_zhpgvx( &itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il,
                       &iu, &abstol, m, w, z, &ldz, work, rwork, iwork, ifail,
                       &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int ncols_z = ( LAPACKE_lsame( range, 'a' ) ||
                             LAPACKE_lsame( range, 'v' ) ) ? n :
                             ( LAPACKE_lsame( range, 'i' ) ? (iu-il+1) : 1);
        lapack_int ldz_t = MAX(1,n);
        lapack_complex_double* z_t = NULL;
        lapack_complex_double* ap_t = NULL;
        lapack_complex_double* bp_t = NULL;
        /* Check leading dimension(s) */
        if( ldz < ncols_z ) {
            info = -17;
            LAPACKE_xerbla( "LAPACKE_zhpgvx_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        if( LAPACKE_lsame( jobz, 'v' ) ) {
            z_t = (lapack_complex_double*)
                LAPACKE_malloc( sizeof(lapack_complex_double) *
                                ldz_t * MAX(1,ncols_z) );
            if( z_t == NULL ) {
                info = LAPACK_TRANSPOSE_MEMORY_ERROR;
                goto exit_level_0;
            }
        }
        ap_t = (lapack_complex_double*)
            LAPACKE_malloc( sizeof(lapack_complex_double) *
                            ( MAX(1,n) * MAX(2,n+1) ) / 2 );
        if( ap_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        bp_t = (lapack_complex_double*)
            LAPACKE_malloc( sizeof(lapack_complex_double) *
                            ( MAX(1,n) * MAX(2,n+1) ) / 2 );
        if( bp_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_2;
        }
        /* Transpose input matrices */
        LAPACKE_zhp_trans( matrix_layout, uplo, n, ap, ap_t );
        LAPACKE_zhp_trans( matrix_layout, uplo, n, bp, bp_t );
        /* Call LAPACK function and adjust info */
        LAPACK_zhpgvx( &itype, &jobz, &range, &uplo, &n, ap_t, bp_t, &vl, &vu,
                       &il, &iu, &abstol, m, w, z_t, &ldz_t, work, rwork, iwork,
                       ifail, &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        if( LAPACKE_lsame( jobz, 'v' ) ) {
            LAPACKE_zge_trans( LAPACK_COL_MAJOR, n, ncols_z, z_t, ldz_t, z,
                               ldz );
        }
        LAPACKE_zhp_trans( LAPACK_COL_MAJOR, uplo, n, ap_t, ap );
        LAPACKE_zhp_trans( LAPACK_COL_MAJOR, uplo, n, bp_t, bp );
        /* Release memory and exit */
        LAPACKE_free( bp_t );
exit_level_2:
        LAPACKE_free( ap_t );
exit_level_1:
        if( LAPACKE_lsame( jobz, 'v' ) ) {
            LAPACKE_free( z_t );
        }
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            LAPACKE_xerbla( "LAPACKE_zhpgvx_work", info );
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_zhpgvx_work", info );
    }
    return info;
}
