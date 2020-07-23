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
#define VARIABLE 1
#define PI 3.14159265358979323846

#if VARIABLE
#define ellipticCoefficient3D(vid,x, y, z, lambda_0, lambda_1)    \
  {                                         \
    if(vid == 0) { \
      lambda_0 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0;            \
      lambda_1 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0;            \
    } \
    if(vid == 1) { \
      lambda_0 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0;            \
      lambda_1 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0;            \
    } \
    if(vid == 2) { \
      lambda_0 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0;            \
      lambda_1 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0;            \
    } \
  }

/* forcing function   */
#define ellipticForcing3D(vid,x, y, z, lambda, f)  \
  {                                         \
    if(vid == 0) { \
      dfloat sxy  = sin(PI * x) * sin(PI * y); \
      dfloat sxz  = sin(PI * x) * sin(PI * z); \
      dfloat syz  = sin(PI * y) * sin(PI * z); \
      dfloat sxyz  = sin(PI * x) * sin(PI * y) * sin(PI * z); \
      f = sxyz * (1.0 + 3.0 * PI * PI + 1.0 * sxyz + 6.0 * PI * PI * sxyz) - PI * PI * \
          (sxy * sxy + sxz * sxz + syz * syz); \
    } \
    if(vid == 1) { \
      dfloat sxy  = sin(PI * x) * sin(PI * y); \
      dfloat sxz  = sin(PI * x) * sin(PI * z); \
      dfloat syz  = sin(PI * y) * sin(PI * z); \
      dfloat sxyz  = sin(PI * x) * sin(PI * y) * sin(PI * z); \
      f = sxyz * (1.0 + 3.0 * PI * PI + 1.0 * sxyz + 6.0 * PI * PI * sxyz) - PI * PI * \
          (sxy * sxy + sxz * sxz + syz * syz); \
    } \
    if(vid == 2) { \
      dfloat sxy  = sin(PI * x) * sin(PI * y); \
      dfloat sxz  = sin(PI * x) * sin(PI * z); \
      dfloat syz  = sin(PI * y) * sin(PI * z); \
      dfloat sxyz  = sin(PI * x) * sin(PI * y) * sin(PI * z); \
      f = sxyz * (1.0 + 3.0 * PI * PI + 1.0 * sxyz + 6.0 * PI * PI * sxyz) - PI * PI * \
          (sxy * sxy + sxz * sxz + syz * syz); \
    } \
  }

/* Dirichlet boundary condition   */
#define ellipticDirichletCondition3D(t,vid,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    if(vid == 0) {   \
      uB  = sin(PI * x) * sin(PI * y) * sin(PI * z);   \
      uxB = uxM;   \
      uyB = uyM;   \
      uzB = uzM;   \
    } \
    if(vid == 1) { \
      uB  = sin(PI * x) * sin(PI * y) * sin(PI * z);   \
      uxB = uxM;   \
      uyB = uyM;   \
      uzB = uzM;   \
    } \
    if(vid == 2) { \
      uB  = sin(PI * x) * sin(PI * y) * sin(PI * z);   \
      uxB = uxM;   \
      uyB = uyM;   \
      uzB = uzM;   \
    } \
  }

/* Neumann boundary condition   */
#define ellipticNeumannCondition3D(t, vid, x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    if(vid == 0) {   \
      uB  = uM;    \
      dfloat lambda_0 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0; \
      uxB = -lambda_0 * PI * cos(PI * x) * sin(PI * y) * sin(PI * z);   \
      uyB = -lambda_0 * PI * sin(PI * x) * cos(PI * y) * sin(PI * z);   \
      uzB = -lambda_0 * PI * sin(PI * x) * sin(PI * y) * cos(PI * z);   \
    } \
    if(vid == 1) {   \
      uB  = uM;    \
      dfloat lambda_0 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0; \
      uxB = -lambda_0 * PI * cos(PI * x) * sin(PI * y) * sin(PI * z);   \
      uyB = -lambda_0 * PI * sin(PI * x) * cos(PI * y) * sin(PI * z);   \
      uzB = -lambda_0 * PI * sin(PI * x) * sin(PI * y) * cos(PI * z);   \
    } \
    if(vid == 2) {   \
      uB  = uM;    \
      dfloat lambda_0 = sin(PI * x) * sin(PI * y) * sin(PI * z) + 1.0; \
      uxB = -lambda_0 * PI * cos(PI * x) * sin(PI * y) * sin(PI * z);   \
      uyB = -lambda_0 * PI * sin(PI * x) * cos(PI * y) * sin(PI * z);   \
      uzB = -lambda_0 * PI * sin(PI * x) * sin(PI * y) * cos(PI * z);   \
    } \
  }
#else // if VARIABLE
/* forcing function   */
#define ellipticForcing3D(vid,x, y, z, lambda, f)  \
  {                                         \
    if(vid == 0) { \
      f  = (3 * PI * PI + lambda) * sin(PI * x) * sin(PI * y) * sin(PI * z);   \
    } \
    if(vid == 1) { \
      f  = (3 * PI * PI + lambda) * sin(PI * x) * sin(PI * y) * sin(PI * z);   \
    } \
    if(vid == 2) { \
      f  = (3 * PI * PI + lambda) * sin(PI * x) * sin(PI * y) * sin(PI * z);   \
    } \
  }

/* Dirichlet boundary condition   */
#define ellipticDirichletCondition3D(t,vid,x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    if(vid == 0) {   \
      uB  = sin(PI * x) * sin(PI * y) * sin(PI * z);   \
      uxB = uxM;   \
      uyB = uyM;   \
      uzB = uzM;   \
    } \
    if(vid == 1) { \
      uB  = sin(PI * x) * sin(PI * y) * sin(PI * z);   \
      uxB = uxM;   \
      uyB = uyM;   \
      uzB = uzM;   \
    } \
    if(vid == 2) { \
      uB  = sin(PI * x) * sin(PI * y) * sin(PI * z);   \
      uxB = uxM;   \
      uyB = uyM;   \
      uzB = uzM;   \
    } \
  }

/* Neumann boundary condition   */
#define ellipticNeumannCondition3D(t, vid, x,y,z,nx,ny,nz,uM,uxM,uyM,uzM,uB,uxB,uyB,uzB)  \
  {              \
    if(vid == 0) {   \
      uB  = uM;    \
      uxB = -PI * cos(PI * x) * sin(PI * y) * sin(PI * z);   \
      uyB = -PI * sin(PI * x) * cos(PI * y) * sin(PI * z);   \
      uzB = -PI * sin(PI * x) * sin(PI * y) * cos(PI * z);   \
    } \
    if(vid == 1) {   \
      uB  = uM;    \
      uxB = -PI * cos(PI * x) * sin(PI * y) * sin(PI * z);   \
      uyB = -PI * sin(PI * x) * cos(PI * y) * sin(PI * z);   \
      uzB = -PI * sin(PI * x) * sin(PI * y) * cos(PI * z);   \
    } \
    if(vid == 2) {   \
      uB  = uM;    \
      uxB = -PI * cos(PI * x) * sin(PI * y) * sin(PI * z);   \
      uyB = -PI * sin(PI * x) * cos(PI * y) * sin(PI * z);   \
      uzB = -PI * sin(PI * x) * sin(PI * y) * cos(PI * z);   \
    } \
  }

#endif
