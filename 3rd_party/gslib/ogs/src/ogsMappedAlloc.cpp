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

#include "ogstypes.h"
#include "ogs.hpp"

void *ogsHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem, occa::memory &h_mem){

#if 0

  mem = device.malloc(size, source);
  
  h_mem = device.mappedAlloc(size, source);
  
  void *ptr = h_mem.getMappedPointer();
#endif
  
  occa::properties props;
  props["mapped"] = true;
  
  if(source!=NULL)
    mem =  device.malloc(size, source);
  else
    mem =  device.malloc(size);

  h_mem =  device.malloc(size, props);
  
  void *ptr = h_mem.ptr(props);
  
  
 return ptr;


  

}
