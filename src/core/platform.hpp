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
#ifndef _platform_hpp_
#define _platform_hpp_
#include "nrssys.hpp"
class platform_t final {
private:
  occa::device device;
  occa::properties kernelInfo;
  MPI_Comm comm;
  static platform_t* singleton;
  platform_t(occa::device& _device, MPI_Comm& _comm)
  : device(_device),
  comm(_comm)
  {
    setup();
  }
  ~platform_t(){
  }
  void setup();
public:
  occa::device& getDevice(){ return device; }
  occa::properties& getKernelInfo(){ return kernelInfo; }
  static platform_t * initialize(occa::device& _device, MPI_Comm& _comm){
    if(!singleton){
      singleton = new platform_t(_device, _comm);
    }
    return singleton;
  }
  static platform_t * getSingleton(){
    return singleton;
  }
};

#endif