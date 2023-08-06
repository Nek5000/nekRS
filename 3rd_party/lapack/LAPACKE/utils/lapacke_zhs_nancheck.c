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

/* Check a matrix for NaN entries. */

lapack_logical LAPACKE_zhs_nancheck( int matrix_layout, lapack_int n,
                                      const lapack_complex_double *a,
                                      lapack_int lda )
{
    lapack_logical subdiag_nans;

    if( a == NULL ) return (lapack_logical) 0;

    /* Check subdiagonal first */
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        subdiag_nans = LAPACKE_z_nancheck( n-1, &a[1], lda+1 );
    } else if ( matrix_layout == LAPACK_ROW_MAJOR ) {
        subdiag_nans = LAPACKE_z_nancheck( n-1, &a[lda], lda+1 );
    } else {
        return (lapack_logical) 0;
    }

    /* Check upper triangular if subdiagonal has no NaNs. */
    return subdiag_nans || LAPACKE_ztr_nancheck( matrix_layout, 'u', 'n',
                                                 n, a, lda);
}
