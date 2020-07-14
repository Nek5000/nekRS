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

#ifndef OCCA_TIMER_HEADER
#define OCCA_TIMER_HEADER

#include "occa.hpp"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <stack>
#include <map>
#include <iomanip>
#include <utility>
#include <algorithm>

namespace occa {

  double currentTime();
  
  class timerTraits{
  public:
    double timeTaken;
    double selfTime;
    int    numCalls;
    double flopCount;
    double bandWidthCount;
    int treeDepth;
    std::vector<std::string> childs;

    timerTraits();
  };

  class timer{

    bool profileKernels;
    bool profileApplication;
    bool deviceInitialized;

    occa::device occaHandle;

  public:

    timer();

    void initTimer(const occa::device &deviceHandle);

    // NBN: allow toggle from menu
    inline void setKernelProfiling(bool b) { profileKernels = b; }
    inline void setApplicationProfiling(bool b) { profileApplication = b; }

    std::stack<std::string> keyStack;
    std::stack<double> timeStack;

    std::map<std::stack<std::string>, timerTraits> times;

    void tic(std::string key);

    double toc(std::string key);

    double toc(std::string key, double flops);

    double toc(std::string key, occa::kernel &kernel);

    double toc(std::string key, occa::kernel &kernel, double flops);

    double toc(std::string key, double flops, double bw);

    double toc(std::string key, occa::kernel &kernel, double flops, double bw);

    double print_recursively(std::vector<std::string> &childs,
                             double parentTime,
                             double overallTime);

    // struct myclass {
    //   bool operator() (std::pair<std::string, timerTraits> &a,
    // 		     std::pair<std::string, timerTraits> &b){
    //     return(a.second.selfTime > b.second.selfTime);
    //   }
    // } compareSelfTimes;


    void printTimer();
  };


  extern timer globalTimer;

  extern double dataTransferred;

  void initTimer(const occa::device &deviceHandle);

  void tic(std::string key);

  double toc(std::string key);

  double toc(std::string key, occa::kernel &kernel);

  double toc(std::string key, double fp);

  double toc(std::string key, occa::kernel &kernel, double fp);

  double toc(std::string key, double fp, double bw);

  double toc(std::string key, occa::kernel &kernel, double fp, double bw);

  void printTimer();
}

extern "C" { // Start C Linkage
void occaTimerTic(occa::device device,std::string name);
void occaTimerToc(occa::device device,std::string name);
} // End C Linkage


#endif

