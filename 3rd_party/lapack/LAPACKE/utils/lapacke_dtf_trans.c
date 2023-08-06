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
* Contents: Native C interface to LAPACK utility function
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

/* Converts input RFP matrix from row-major(C) to column-major(Fortran)
 * layout or vice versa.
 * This functions does copy diagonal for both unit and non-unit cases.
 */

void LAPACKE_dtf_trans( int matrix_layout, char transr, char uplo, char diag,
                        lapack_int n, const double *in,
                        double *out )
{
    lapack_int row, col;
    lapack_logical rowmaj, ntr, lower, unit;

    if( in == NULL || out == NULL ) return ;

    rowmaj = (matrix_layout == LAPACK_ROW_MAJOR);
    ntr    = LAPACKE_lsame( transr, 'n' );
    lower  = LAPACKE_lsame( uplo,   'l' );
    unit   = LAPACKE_lsame( diag,   'u' );

    if( ( !rowmaj && ( matrix_layout != LAPACK_COL_MAJOR ) ) ||
        ( !ntr    && !LAPACKE_lsame( transr, 't' ) &&
                     !LAPACKE_lsame( transr, 'c' ) ) ||
        ( !lower  && !LAPACKE_lsame( uplo,   'u' ) ) ||
        ( !unit   && !LAPACKE_lsame( diag,   'n' ) ) ) {
        /* Just exit if input parameters are wrong */
        return;
    }

    /* Determine parameters of array representing RFP */
    if( ntr ) {
        if( n%2 == 0 ) {
            row = n + 1;
            col = n / 2;
        } else {
            row = n;
            col = (n + 1) / 2;
        }
    } else {
        if( n%2 == 0 ) {
            row = n / 2;
            col = n + 1;
        } else {
            row = (n + 1) / 2;
            col = n;
        }
    }

    /* Perform conversion: */
    if( rowmaj ) {
        LAPACKE_dge_trans( LAPACK_ROW_MAJOR, row, col, in, col, out, row );
    } else {
        LAPACKE_dge_trans( LAPACK_COL_MAJOR, row, col, in, row, out, col );
    }
}
