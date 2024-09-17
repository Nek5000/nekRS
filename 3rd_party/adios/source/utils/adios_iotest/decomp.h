/*
 * decomp.h
 *
 *  Created on: Oct 2018
 *      Author: Norbert Podhorszki
 */

#ifndef DECOMP_H_
#define DECOMP_H_

#include <stddef.h> // size_t

/** Decompose M processes with N-dimensional decomposition in a Column major
 * fashion.
 * INPUT: ndim as N dimension
 *        rank as the MPI rank of the calling process
 *        decomp[] array for the number of processes in each dimension
 * OUTPUT: pos[] array describes the position of this rank in each dimension
 *
 * Requirement: pos[] must be pre-allocated to hold ndim elements
 *
 * Examples:
 * 1. 2-dim, 6 processes with 2x3 decomposition
 *     ndim = 2, decomp[]={2,3}
 *     decomposition calculated according to this table
 *                    pos[1]
 *                  0   1   2
 *                 ----------
 *     pos[0]  0 |  0   2   4
 *             1 |  1   3   5
 *
 */

void decompColumnMajor(const size_t ndim, const size_t rank, const size_t *decomp, size_t *pos);

/** Decompose M processes with N-dimensional decomposition in a Row major
 * fashion.
 * INPUT: ndim as N dimension
 *        rank as the MPI rank of the calling process
 *        decomp[] array for the number of processes in each dimension
 * OUTPUT: pos[] array describes the position of this rank in each dimension
 *
 * Requirement: pos[] must be pre-allocated to hold ndim elements
 *
 * Examples:
 * 1. 2-dim, 6 processes with 2x3 decomposition
 *     ndim = 2, decomp[]={2,3}
 *     decomposition calculated according to this table
 *                    pos[1]
 *                  0   1   2
 *                 ----------
 *     pos[0]  0 |  0   1   2
 *             1 |  3   4   5
 *
 */

void decompRowMajor(const size_t ndim, const size_t rank, const size_t *decomp, size_t *pos);

#endif /* DECOMP_H */
