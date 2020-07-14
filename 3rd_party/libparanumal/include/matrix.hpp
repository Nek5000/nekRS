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

#ifndef __MATRIX
#define __MATRIX

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <execinfo.h>

static int clarrayrequest = 0;

template <class T>
class matrix {

private:
  int Nrows;
  int Ncolumns;
  T *data;

  void allocate();

public:

  matrix();

  matrix(int nr, int nc);

  // copy constructor
  matrix(const matrix <T> &A);

  // assignment operator
  matrix <T> operator= (const matrix <T> &A);

  void operator += (const matrix <T> &A);

  // assignment operator from double **array (assumes initialized)
  matrix <T> operator= (const T *A);

  // assignment operator (all values set to d)
  matrix <T> operator= (const double &d);

  ~matrix();

  int nrows() const;

  int ncolumns() const;

  T *c_array();

  void resize(int nr, int nc);

  void resize(int nr);

  matrix <T> transpose();

  // column major serial access - 1-indexed
  T operator[] (int r) const ;

  // column major serial access - 1-indexed
  T & operator[] (int r) ;

  T operator() (int r, int c) const ;

  T & operator() (int r, int c) ;

  T operator() (int r) const;

  T & operator() (int r) ;

  // permute
  matrix <T> operator[] (const matrix <int> &ind)  ;

  // in-place sort in ascending order
  void slow_sort(int index);

  // member function to randomize entries
  void randomize();

  // sort columns using user supplied comparison function
  void sort( int (*compare)(const void *, const void *) );

  void symeig(matrix <T> &d, matrix <T> &v);

  void eig(matrix <T> &WR, matrix <T> &WI, matrix <T> &VL, matrix <T> &VR);

  int byteCount() const;

  int entryCount() const;

  int size();

  int size() const;

  long double frobenius() const;

  T maxentry() const;

  T minentry() const;

  matrix <T> inverse();

};

using namespace std;

template <class T>
ostream & operator << (ostream &os, matrix <T> & A);

template <class T>
matrix <T> operator* (const matrix <T> & A, const matrix <T> &B);

template <class T>
matrix <T> operator+ (const matrix <T> & A, const matrix <T> &B);

template <class T>
matrix <T> operator- (const matrix <T> & A, const matrix <T> &B);

// general left matrix inverse not implemented
template <class T>
matrix <T> operator| (const matrix <T> & A, const matrix <T> &B);

// use #define to create a "compiler substitution rule"
#define imatrix matrix<int>
#define dmatrix matrix<double>


#include "matrix.tpp"

#endif
