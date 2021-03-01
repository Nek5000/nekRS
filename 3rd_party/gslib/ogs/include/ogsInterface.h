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

#ifndef OGS_INTERFACE_H
#define OGS_INTERFACE_H 1

extern "C"
{
  void *ogsHostSetup(MPI_Comm comm, dlong Ngather, hlong *gatherIds, int unique, int verbose);
  void  ogsGsUnique(hlong *gatherIds, dlong Ngather, MPI_Comm comm);

  void ogsHostGatherScatter    (void *v, const char *type, const char *op, void *gsh);
  void ogsHostGatherScatterVec (void *v, const int k, const char *type, const char *op, void *gsh);
  void ogsHostGatherScatterMany(void *v, const int k, const char *type, const char *op, void *gsh);
  
  void ogsHostGather    (void *v, const char *type, const char *op, void *gsh);
  void ogsHostGatherVec (void *v, const int k, const char *type, const char *op, void *gsh);
  void ogsHostGatherMany(void *v, const int k, const char *type, const char *op, void *gsh);
  
  void ogsHostScatter    (void *v, const char *type, const char *op, void *gsh);
  void ogsHostScatterVec (void *v, const int k, const char *type, const char *op, void *gsh);
  void ogsHostScatterMany(void *v, const int k, const char *type, const char *op, void *gsh);
  
  void ogsHostFree(void *gsh);

}

#endif