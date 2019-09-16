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

#define WALL 1
#define INLET 2
#define OUTLET 3
#define XSLIP 4
#define YSLIP 5
#define ZSLIP 6

struct bcData
{
   int idM;
   int fieldOffset;
   int id;

   dfloat time;
   dfloat x, y, z;
   dfloat nx, ny, nz;

   dfloat uM, vM, wM;
   dfloat uP, vP, wP;
   dfloat uxM, uyM, uzM;
   dfloat vxM, vyM, vzM;
   dfloat wxM, wyM, wzM;
   dfloat uxP, uyP, uzP;
   dfloat vxP, vyP, vzP;
   dfloat wxP, wyP, wzP;

   dfloat pM;
   dfloat pP, pxP, pyP, pzP;

   dfloat* wrk;

};
