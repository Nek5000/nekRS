/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Created by Dmitry Ganyushin ganyushindi@ornl.gov
 */

#ifndef ADIOS2EXAMPLES_DECOMP_H
#define ADIOS2EXAMPLES_DECOMP_H

#include <stddef.h>

extern size_t get_random(int minv, int maxv, int rank);
extern void gather_decomp_1d(size_t *, size_t *, size_t *);
extern void decomp_1d(size_t *, size_t *, size_t *);
#endif // ADIOS2EXAMPLES_DECOMP_H
