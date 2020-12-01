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

using std::stringstream;
using std::fstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;


class setupAide {
private:
  vector<string> data;
  vector<string> keyword;

public:
  setupAide();
  setupAide(string);

  setupAide(const setupAide&);
  setupAide& operator=(const setupAide&);

  string readFile(string);
  void read(string);

  string getArgs(string);

  void setArgs(string key, string value);

  template <class T>
  int getArgs(string, T&);

  template <class T>
  int getArgs(string, vector<T>&);

  int getArgs(string, vector<string>&, string);


  int compareArgs(string key, string token);

  vector<string> &getData(){ return data; }
  vector<string> &getKeyword() { return keyword; }

  friend std::ostream & operator << (std::ostream &out, const setupAide &aide);  
};

#include<setupAide.tpp>


#endif

// Combined header/source needs explicit instantiation
