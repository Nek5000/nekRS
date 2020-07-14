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

#ifndef ELLIPTIC_TIMER_H_
#define ELLIPTIC_TIMER_H_ 1
#include <string.h>
#include "mpi.h"
#include <map>
#include <vector>
#include <utility>
#include <iostream>

//#define ENABLE_TIMER 1
struct Entry final{
  std::string name;
  int level;
  Entry(std::string _name ,int _level) : name(_name), level(_level) {}
  bool operator==(const Entry& other) const
  {
    return (other.name == this->name) && (other.level == this->level);
  }
  bool operator<(const Entry& other) const
  {
    if(this->name < other.name){
        if(this->level < other.level){
            return true;
        }
    }
    return false;
  }
};
class PerformanceTimer final
{
public:
    static PerformanceTimer& getInstance()
    {
        static PerformanceTimer instance;
        return instance;
    }
    PerformanceTimer(PerformanceTimer const&) = delete;
    void operator=(PerformanceTimer const&) = delete;
    void tic(const Entry& e);
    void toc(const Entry& e);
    std::pair<double,double> stats(const Entry& e) const;
    void dump() const;
private:
    PerformanceTimer() {}
    std::map<Entry,std::vector<double>> entry_to_times_map;
};
#ifdef ENABLE_TIMER
#define START(myStrName,myLevel) PerformanceTimer::getInstance().tic(Entry{myStrName,myLevel});
#define STOP(myStrName,myLevel) PerformanceTimer::getInstance().toc(Entry{myStrName,myLevel});
#else
#define START(myStrName,myLevel) do {} while(0);
#define STOP(myStrName,myLevel) do {} while(0);
#endif
#endif

