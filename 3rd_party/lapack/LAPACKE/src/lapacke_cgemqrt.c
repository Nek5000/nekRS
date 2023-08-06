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
******************************************************************************
* Contents: Native high-level C interface to LAPACK function cgemqrt
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_cgemqrt( int matrix_layout, char side, char trans,
                            lapack_int m, lapack_int n, lapack_int k,
                            lapack_int nb, const lapack_complex_float* v,
                            lapack_int ldv, const lapack_complex_float* t,
                            lapack_int ldt, lapack_complex_float* c,
                            lapack_int ldc )
{
    lapack_int nrows_v;
    lapack_int info = 0;
    lapack_complex_float* work = NULL;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        LAPACKE_xerbla( "LAPACKE_cgemqrt", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        nrows_v = LAPACKE_lsame( side, 'L' ) ? m :
                             ( LAPACKE_lsame( side, 'R' ) ? n : 0 );
        if( LAPACKE_cge_nancheck( matrix_layout, m, n, c, ldc ) ) {
            return -12;
        }
        if( LAPACKE_cge_nancheck( matrix_layout, nb, k, t, ldt ) ) {
            return -10;
        }
        if( LAPACKE_cge_nancheck( matrix_layout, nrows_v, k, v, ldv ) ) {
            return -8;
        }
    }
#endif
    /* Allocate memory for working array(s) */
    work = (lapack_complex_float*)
        LAPACKE_malloc( sizeof(lapack_complex_float) * MAX(1,m) * MAX(1,nb) );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    /* Call middle-level interface */
    info = LAPACKE_cgemqrt_work( matrix_layout, side, trans, m, n, k, nb, v, ldv,
                                 t, ldt, c, ldc, work );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        LAPACKE_xerbla( "LAPACKE_cgemqrt", info );
    }
    return info;
}
