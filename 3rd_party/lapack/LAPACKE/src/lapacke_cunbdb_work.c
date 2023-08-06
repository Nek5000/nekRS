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
* Contents: Native middle-level C interface to LAPACK function cunbdb
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int LAPACKE_cunbdb_work( int matrix_layout, char trans, char signs,
                                lapack_int m, lapack_int p, lapack_int q,
                                lapack_complex_float* x11, lapack_int ldx11,
                                lapack_complex_float* x12, lapack_int ldx12,
                                lapack_complex_float* x21, lapack_int ldx21,
                                lapack_complex_float* x22, lapack_int ldx22,
                                float* theta, float* phi,
                                lapack_complex_float* taup1,
                                lapack_complex_float* taup2,
                                lapack_complex_float* tauq1,
                                lapack_complex_float* tauq2,
                                lapack_complex_float* work, lapack_int lwork )
{
    lapack_int info = 0;
    /* LAPACK function works with matrices in both layouts. It is supported
     * through TRANS parameter. So all conversion between layouts can be
     * completed in LAPACK function. See the table below which describes how
     * every LAPACKE call is forwarded to corresponding LAPACK call.
     *
     *   matrix_layout    | trans_LAPACKE | -> trans_LAPACK
     *                    |   (trans)     |       (ltrans)
     *   -----------------+---------------+----------------
     *   LAPACK_COL_MAJOR |     'N'       | ->     'N'
     *   LAPACK_COL_MAJOR |     'T'       | ->     'T'
     *   LAPACK_ROW_MAJOR |     'N'       | ->     'T'
     *   LAPACK_ROW_MAJOR |     'T'       | ->     'T'
     *       (note that for row major layout trans parameter is ignored)
     */
    if( matrix_layout == LAPACK_COL_MAJOR ||
        matrix_layout == LAPACK_ROW_MAJOR ) {
        char ltrans;
        if( !LAPACKE_lsame( trans, 't' ) && matrix_layout == LAPACK_COL_MAJOR ) {
            ltrans = 'n';
        } else {
            ltrans = 't';
        }
        /* Call LAPACK function and adjust info */
        LAPACK_cunbdb( &ltrans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12,
                       x21, &ldx21, x22, &ldx22, theta, phi, taup1, taup2,
                       tauq1, tauq2, work, &lwork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else {
        info = -1;
        LAPACKE_xerbla( "LAPACKE_cunbdb_work", info );
    }
    return info;
}
