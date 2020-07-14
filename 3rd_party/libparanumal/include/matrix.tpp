#include<assert.h>

#include "trace.hpp"
#include "matrix.hpp"

template<class T>
matrix<T>::matrix(){

  Nrows = 0;
  Ncolumns = 0;
  data = 0;

}

template<class T>
matrix<T>::matrix(int nr, int nc){

  Nrows = nr;
  Ncolumns = nc;

  allocate();

}

// copy constructor
template<class T>
matrix<T>:: matrix(const matrix <T> &A){

  Nrows = A.Nrows;
  Ncolumns = A.Ncolumns;

  allocate();

  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = A(r,c);

}

template<class T>
void matrix<T>::allocate(){
  if(Nrows*Ncolumns)
    data = new T[Nrows*Ncolumns]();
  else {
    Nrows = 0;
    Ncolumns = 0;
    data = NULL;
  }
}

// assignment operator
template<class T> matrix<T>
matrix<T>::operator= (const matrix <T> &A){

  if(this != &A){
    if(Nrows != A.Nrows || Ncolumns != A.Ncolumns){

      if(data)
	delete [] data;

      Nrows = A.Nrows;
      Ncolumns = A.Ncolumns;

      allocate();

    }

    for(int r=1; r<=Nrows; r=r+1)
      for(int c=1; c<=Ncolumns; c=c+1)
	(*this)(r,c) = A(r,c);

  }
  return *this;

}

template<class T>
void matrix<T>::operator += (const matrix <T> &A){

  if(this != &A){
    if(Nrows != A.Nrows || Ncolumns != A.Ncolumns){

      if(data)
	delete [] data;

      Nrows = A.Nrows;
      Ncolumns = A.Ncolumns;

      allocate();

    }

    for(int r=1; r<=Nrows; r=r+1)
      for(int c=1; c<=Ncolumns; c=c+1)
	(*this)(r,c) += A(r,c);

  }
}


// assignment operator from double **array (assumes initialized)
template<class T>  matrix <T>
matrix<T>::operator= (const T *A){

  // load from C style matrix
  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = A[Ncolumns*(r-1) + c-1];

  return *this;

}

// assignment operator (all values set to d)
template<class T> matrix <T>
matrix<T>::operator= (const double &d){

  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = d;

  return *this;
}

template<class T>
matrix<T>:: ~matrix(){
  if(data)
    delete [] data;
}

template<class T> int
matrix<T>:: nrows() const{
  return Nrows;
}

template<class T> int
matrix<T>::ncolumns() const{
  return Ncolumns;
}

template<class T>
T* matrix<T>::c_array(){
  return data;
}

template<class T> void
matrix<T>::resize(int nr){
  resize(nr,1);
}


template<class T> void
matrix<T>::resize(int nr, int nc){

  // temporary
  matrix <T> dupe(nr,nc);

  dupe = 0.;

  if(data){
    // entry by entry copy to dupe
    for(int r=1; r<=nr; r=r+1)
      for(int c=1; c<=nc; c=c+1)
	if(r<=Nrows && c<=Ncolumns)
	  dupe(r,c) = (*this)(r,c);

    delete [] data;
  }

  Nrows = nr;
  Ncolumns = nc;

  allocate();

  for(int r=1; r<=nr; r=r+1)
    for(int c=1; c<=nc; c=c+1)
      (*this)(r,c) = dupe(r,c);
}

template<class T> matrix <T>
matrix<T>::transpose(){

  // temporary
  matrix <T> AT(Ncolumns, Nrows);

  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      AT(c,r) = (*this)(r,c);

  return AT;
}

  // column major serial access - 1-indexed
template<class T>
T matrix<T>::operator[] (int r) const {

  // bounds check here [ could also throw an exception ]

  assert(r>0);
  assert(r<=Nrows*Ncolumns);
  if(r<=0 || r> Nrows*Ncolumns){
    std::cout << "matrix::operator[] Out of bounds access" << std::endl;
    trace(std::cout, 20);
    exit(-1);
  }

  return data[ (r-1) ];
}

  // column major serial access - 1-indexed
template<class T>
T & matrix<T>::operator[] (int r) {

  assert(r>0);
  assert(r<=Nrows*Ncolumns);
  // bounds check here [ could also throw an exception ]
  if(r<=0 || r> Nrows*Ncolumns){
    std::cout << "matrix::operator[] Out of bounds access" << std::endl;
    trace(std::cout, 20);
    exit(-1);
  }

  return data[ (r-1) ];
}

template<class T>
T matrix<T>::operator() (int r, int c) const {

  assert(r>0);
  assert(r<=Nrows);

  assert(c>0);
  assert(c<=Ncolumns);

  // bounds check here [ could also throw an exception ]
  if(r<=0 || r>Nrows || c<=0 || c>Ncolumns){
    std::cout << "matrix::operator() Out of bounds [" << nrows() << "," << ncolumns()
	      << "] access" << "(r,c)=(" << r << "," << c << ")" << std::endl;
    trace(std::cout, 20);
    //      exit(-1);
  }

  return data[ (r-1) + (c-1)*Nrows ];
}

template<class T>
T & matrix<T>::operator() (int r, int c) {

  assert(r>0);
  assert(r<=Nrows);

  assert(c>0);
  assert(c<=Ncolumns);
  // bounds check here [ could also throw an exception ]
  if(r<=0 || r>Nrows || c<=0 || c>Ncolumns){
    std::cout << "matrix::operator() Out of bounds [" << nrows() << "," << ncolumns()
	      << "] access" << "(r,c)=(" << r << "," << c << ")" << std::endl;
    trace(std::cout, 20);
    //      exit(-1);
  }

  return data[ (r-1) + (c-1)*Nrows ];
}

template<class T>
T matrix<T>::operator() (int r) const{
  assert(r>0);
  assert(r<=Nrows*Ncolumns);
  // bounds check here [ could also throw an exception ]
  if(r<=0 || r>Nrows*Ncolumns){
    std::cout << "matrix::operator() Out of bounds access" << std::endl;
    trace(std::cout, 20);
    //      exit(-1);
  }

  return data[ r-1 ];
}

template<class T>
T & matrix<T>::operator() (int r) {
  assert(r>0);
  assert(r<=Nrows*Ncolumns);
  // bounds check here [ could also throw an exception ]
  if(r<=0 || r>Nrows*Ncolumns){
    std::cout << "matrix::operator() Out of bounds access" << std::endl;
    trace(std::cout, 20);
    //      exit(-1);
  }

  return data[ r-1 ];
}

  // permute
template<class T> matrix <T>
matrix<T>::operator[] (const matrix <int> &ind)  {

  matrix <T> C(ind.nrows(), ind.ncolumns());

  for(int r=1;r<=ind.nrows();r=r+1){
    for(int c=1;c<=ind.ncolumns();c=c+1){
      C(r,c) = (*this)[ind(r,c)];
    }
  }
  return C;
}


  // in-place sort in ascending order
template<class T>
void matrix<T>::slow_sort(int index){

  // O( (4*K)^2 )

  // slow sort
  for(int r1=1;r1<=Nrows;r1=r1+1){
    for(int r2=r1+1;r2<=Nrows;r2=r2+1){
      if( (*this)(r2,index) < (*this)(r1,index) ){
	for(int c=1;c<=Ncolumns;c=c+1){
	  // switch rows r1 and r2
	  T t = (*this)(r1,c);
	  (*this)(r1,c) = (*this)(r2,c);
	  (*this)(r2,c) = t;
	}
      }
    }
  }
}

// member function to randomize entries
template<class T>
void matrix<T>::randomize(){
  int r, c;
  for(c=1;c<=Ncolumns;++c){
    for(r=1;r<=Nrows;++r){
      (*this)(r,c) = drand48();
    }
  }
}


// sort columns using user supplied comparison function
template<class T>
void matrix<T>::sort( int (*compare)(const void *, const void *) ){

  qsort(data, Ncolumns, Nrows*sizeof(T), compare);

}

// template<class T>
// void matrix<T>::symeig(matrix <T> &d, matrix <T> &v){

//   // not implemented for general T
// }


// template<class T>
// void matrix<T>::eig(matrix <T> &WR, matrix <T> &WI, matrix <T> &VL, matrix <T> &VR){

//   // not implemented for general T
//   // beware assume eigenvecs are real
// }

template<class T>
int matrix<T>::byteCount() const{

  return sizeof(T)*Nrows*Ncolumns;

}

template<class T>
int matrix<T>::entryCount() const{

  return Nrows*Ncolumns;

}


template<class T>
int matrix<T>::size(){
  return Nrows*Ncolumns;
}

template<class T>
int matrix<T>::size() const{
  return Nrows*Ncolumns;
}

template<class T>
long double matrix<T>::frobenius() const{

  long double f(0);

  for(int c=1;c<=Ncolumns;++c)
    for(int r=1;r<=Nrows;++r)
      f += (*this)(r,c)*(*this)(r,c);

  return f;
}

template<class T>
T matrix<T>::maxentry() const{

  T maxval = (*this)(1,1);

  for(int c=1;c<=Ncolumns;++c){
    for(int r=1;r<=Nrows;++r){
      T val = (*this)(r,c);
      maxval = mymax(maxval, val);
    }
  }

  return maxval;
}

template<class T>
T matrix<T>::minentry() const{

  T minval = (*this)(1,1);

  for(int c=1;c<=Ncolumns;++c){
    for(int r=1;r<=Nrows;++r){
      T val = (*this)(r,c);
      minval = mymin(minval, val);
    }
  }

  return minval;
}

template<class T> matrix <T>
matrix<T>::inverse(){

  matrix <T> I(Nrows, Ncolumns);
  for(int n=1;n<=Nrows;++n)
    I(n,n) = 1;

  matrix <T> invA = operator| (*this, I);
  return invA;
}




template <class T>
ostream & operator << (ostream &os, matrix <T> & A){
  int r,c;

  os << std::scientific << "[";
  for(r=1;r<=A.nrows();r=r+1){

    os << "[";

    for(c=1;c<A.ncolumns();c=c+1)
      os << A(r,c) << ", \t";

    os << A(r,c) << "]";

    if(r<A.nrows())
      os << std::endl;
  }
  os << "];" << std::endl;

  return os;

}


template <class T>
matrix <T> operator* (const matrix <T> & A, const matrix <T> &B){

  matrix <T> C(A.nrows(), B.ncolumns());
  //  assert(A.ncolumns() == B.nrows());
  for(int r=1;r<=A.nrows();r=r+1){
    for(int c=1;c<=B.ncolumns();c=c+1){
      // carefully use this to start row-column product
      T s = A(r,1)*B(1,c);
      for(int i=2;i<=A.ncolumns();i=i+1)
	s = s + A(r,i)*B(i,c);
      C(r,c) = s;
    }
  }

  return C;

}


template <class T>
matrix <T> operator+ (const matrix <T> & A, const matrix <T> &B){

  matrix <T> C = A;

  for(int r=1;r<=A.nrows();r=r+1)
    for(int c=1;c<=B.ncolumns();c=c+1)
      C(r,c) = C(r,c)+B(r,c);

  return C;

}

template <class T>
matrix <T> operator- (const matrix <T> & A, const matrix <T> &B){

  matrix <T> C = A;

  for(int r=1;r<=A.nrows();r=r+1)
    for(int c=1;c<=B.ncolumns();c=c+1)
      C(r,c) = C(r,c)-B(r,c);

  return C;

}


