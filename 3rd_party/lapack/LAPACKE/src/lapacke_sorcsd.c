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
* Contents: Native high-level C interface to LAPACK function sorcsd
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_sorcsd( int matrix_layout, char jobu1, char jobu2,
                           char jobv1t, char jobv2t, char trans, char signs,
                           lapack_int m, lapack_int p, lapack_int q, float* x11,
                           lapack_int ldx11, float* x12, lapack_int ldx12,
                           float* x21, lapack_int ldx21, float* x22,
                           lapack_int ldx22, float* theta, float* u1,
                           lapack_int ldu1, float* u2, lapack_int ldu2,
                           float* v1t, lapack_int ldv1t, float* v2t,
                           lapack_int ldv2t )
{
    lapack_int info = 0;
    lapack_int lwork = -1;
    lapack_int* iwork = NULL;
    float* work = NULL;
    float work_query;
    int lapack_layout;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        LAPACKE_xerbla( "LAPACKE_sorcsd", -1 );
        return -1;
    }
    if( LAPACKE_lsame( trans, 'n' ) && matrix_layout == LAPACK_COL_MAJOR ) {
        lapack_layout = LAPACK_COL_MAJOR;
    } else {
        lapack_layout = LAPACK_ROW_MAJOR;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        if( LAPACKE_sge_nancheck( lapack_layout, p, q, x11, ldx11 ) ) {
            return -11;
        }
        if( LAPACKE_sge_nancheck( lapack_layout, p, m-q, x12, ldx12 ) ) {
            return -13;
        }
        if( LAPACKE_sge_nancheck( lapack_layout, m-p, q, x21, ldx21 ) ) {
            return -15;
        }
        if( LAPACKE_sge_nancheck( lapack_layout, m-p, m-q, x22, ldx22 ) ) {
            return -17;
        }
    }
#endif
    /* Allocate memory for working array(s) */
    iwork = (lapack_int*)LAPACKE_malloc( sizeof(lapack_int) * MAX(1,m-MIN(MIN(p,m-p),MIN(q,m-q))) );
    if( iwork == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    /* Query optimal working array(s) size */
    info = LAPACKE_sorcsd_work( matrix_layout, jobu1, jobu2, jobv1t, jobv2t,
                                trans, signs, m, p, q, x11, ldx11, x12, ldx12,
                                x21, ldx21, x22, ldx22, theta, u1, ldu1, u2,
                                ldu2, v1t, ldv1t, v2t, ldv2t, &work_query,
                                lwork, iwork );
    if( info != 0 ) {
        goto exit_level_1;
    }
    lwork = (lapack_int)work_query;
    /* Allocate memory for work arrays */
    work = (float*)LAPACKE_malloc( sizeof(float) * lwork );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_1;
    }
    /* Call middle-level interface */
    info = LAPACKE_sorcsd_work( matrix_layout, jobu1, jobu2, jobv1t, jobv2t,
                                trans, signs, m, p, q, x11, ldx11, x12, ldx12,
                                x21, ldx21, x22, ldx22, theta, u1, ldu1, u2,
                                ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork,
                                iwork );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_1:
    LAPACKE_free( iwork );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        LAPACKE_xerbla( "LAPACKE_sorcsd", info );
    }
    return info;
}
