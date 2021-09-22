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

#ifndef __SETUPAIDE
#define __SETUPAIDE

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>

#include "nrssys.hpp"

class setupAide {
private:
  std::vector<std::string> data;
  std::vector<std::string> keyword;

public:
  setupAide();
  setupAide(std::string);

  setupAide(const setupAide&);
  setupAide& operator=(const setupAide&);

  std::string readFile(std::string);
  void read(std::string);

  std::string getArgs(std::string);

  void setArgs(std::string key, std::string value);

  template <class T>
  int getArgs(std::string, T&);

  template <class T>
  int getArgs(std::string, std::vector<T>&);

  int getArgs(std::string, std::vector<std::string>&, std::string);


  int compareArgs(std::string key, std::string token);

  std::vector<std::string> &getData(){ return data; }
  std::vector<std::string> &getKeyword() { return keyword; }

  friend std::ostream & operator << (std::ostream &out, const setupAide &aide);  
};

#include<setupAide.tpp>


#endif

// Combined header/source needs explicit instantiation
