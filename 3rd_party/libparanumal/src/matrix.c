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

#include "matrix.hpp"


template<>
void matrix<float>::symeig(matrix <float> &W, matrix <float> &V){

  V = *this;
  W.resize(Nrows,1);

  char JOBZ = 'V';
  char UPLO = 'U';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <float> WORK(LWORK,1);

  int INFO = -999;

  //  ssyev_ (&JOBZ, &UPLO, &N, V.c_array(), &LDA, W.c_array(), WORK.c_array(), &LWORK, &INFO );

  //  cout << "INFO: " << INFO << endl;

}

template<>
void matrix<double>::symeig(matrix <double> &W, matrix <double> &V){

  V = *this;
  W.resize(Nrows,1);

  char JOBZ = 'V';
  char UPLO = 'U';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <double> WORK(LWORK,1);

  int INFO = -999;

  //  dsyev_ (&JOBZ, &UPLO, &N, V.c_array(), &LDA, W.c_array(), WORK.c_array(), &LWORK, &INFO );

  //  cout << "INFO: " << INFO << endl;
}

template<>
void matrix<float>::eig(matrix <float> &WR, matrix <float> &WI,
			matrix <float> &VL, matrix <float> &VR){

  matrix <float> A = *this;
  VL = *this;
  VR = *this;
  WR.resize(Nrows,1);
  WI.resize(Nrows,1);

  char JOBVL = 'V';
  char JOBVR = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <float> WORK(LWORK,1);

  int INFO = -999;

  //  sgeev_ (&JOBVL, &JOBVR, &N, A.c_array(), &LDA, WR.c_array(), WI.c_array(), VL.c_array(), &LDA, VR.c_array(), &LDA, WORK.c_array(), &LWORK, &INFO);


  //  cout << "INFO: " << INFO << endl;

}

template<>
void matrix<double>::eig(matrix <double> &WR, matrix <double> &WI,
			 matrix <double> &VL, matrix <double> &VR){

  matrix <double> A = *this;
  VL = *this;
  VR = *this;
  WR.resize(Nrows,1);
  WI.resize(Nrows,1);

  char JOBVL = 'V';
  char JOBVR = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK = N*N;
  matrix <double> WORK(LWORK,1);

  int INFO = -999;

  //  dgeev_ (&JOBVL, &JOBVR, &N, A.c_array(), &LDA, WR.c_array(), WI.c_array(), VL.c_array(), &LDA, VR.c_array(), &LDA, WORK.c_array(), &LWORK, &INFO);


  //  cout << "INFO: " << INFO << endl;

}


// general left matrix inverse not implemented
template<>
matrix <double> operator| (const matrix <double> & A, const matrix <double> &B){

  matrix <double> C = B;
  matrix <double> Acopy = A;

  int N    = A.nrows();
  int NRHS = B.ncolumns();
  int LDA  = N;
  int LDB  = B.nrows();
  int INFO;

  matrix <int> IPIV(N,1);

  //  cout << "NRHS: " << NRHS << endl;

  /* solve A*X = B for X */
#if 0
  dgesv_( &N,
	  &NRHS,
	  Acopy.c_array(),
	  &LDA,
          IPIV.c_array(),
	  C.c_array(),
	  &LDB,
	  &INFO );
#endif
  
  return C;
}


// general left matrix inverse not implemented
template<>
matrix <float> operator| (const matrix <float> & A, const matrix <float> &B){

  matrix <float> C = B;
  matrix <float> Acopy = A;

  int N    = A.nrows();
  int NRHS = B.ncolumns();
  int LDA  = N;
  int LDB  = B.nrows();
  int INFO;

  matrix <int> IPIV(N,1);

  //  cout << "NRHS: " << NRHS << endl;
  //  cout << "Acopy: " << Acopy << endl;

  /* solve A*X = B for X */
#if 0
  sgesv_( &N,
	  &NRHS,
	  Acopy.c_array(),
	  &LDA,
          IPIV.c_array(),
	  C.c_array(),
	  &LDB,
	  &INFO );
#endif
  if(INFO)
    std::cout << "sgesv: INFO=" << INFO << std::endl;


  return C;
}
