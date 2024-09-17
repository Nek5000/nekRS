/*
 * decomp.h
 *
 *  Created on: Oct 2018
 *      Author: Norbert Podhorszki
 */

#include "decomp.h"

void decompColumnMajor(const size_t ndim, const size_t rank, const size_t *decomp, size_t *pos)
{
    // pos[k] = rank / prod(decomp[i], i=0..k-1) % decomp[k]

    size_t prod = 1;
    for (size_t k = 0; k < ndim; ++k)
    {
        pos[k] = rank / prod % decomp[k];
        prod *= decomp[k];
    }
}

void decompRowMajor(const size_t ndim, const size_t rank, const size_t *decomp, size_t *pos)
{
    // pos[k] = rank / prod(decomp[i], i=k+1..n-1) % decomp[k]

    size_t prod = 1;
    for (int k = static_cast<int>(ndim) - 1; k >= 0; --k)
    {
        pos[k] = rank / prod % decomp[k];
        prod *= decomp[k];
    }
}
