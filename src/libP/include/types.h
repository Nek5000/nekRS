/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


//float data type
#if 0
#define DFLOAT_SINGLE
#define dfloat float
#define MPI_DFLOAT MPI_FLOAT
#define dfloatFormat "%f"
#define dfloatString "float"
#else
#define DFLOAT_DOUBLE
#define dfloat double
#define MPI_DFLOAT MPI_DOUBLE
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif

//smoother float data type
#if 1
#define pfloat float
#define MPI_PFLOAT MPI_FLOAT
#define pfloatFormat "%f"
#define pfloatString "float"
#else
#define pfloat double
#define MPI_PFLOAT MPI_DOUBLE
#define pfloatFormat "%lf"
#define pfloatString "double"
#endif

//host index data type
#if 0
#define hlong int
#define MPI_HLONG MPI_INT
#define hlongFormat "%d"
#define hlongString "int"
#else
#define hlong long long int
#define MPI_HLONG MPI_LONG_LONG_INT
#define hlongFormat "%lld"
#define hlongString "long long int"
#endif

//device index data type
#if 1
#define dlong int
#define MPI_DLONG MPI_INT
#define dlongFormat "%d"
#define dlongString "int"
#else
#define dlong long long int
#define MPI_DLONG MPI_LONG_LONG_INT
#define dlongFormat "%lld"
#define dlongString "long long int"
#endif
