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

#include <math.h>

unsigned int hash(const unsigned int value) {
  
  const int p[8] = {
    102679, 102701, 102761, 102763,
    102769, 102793, 102797, 102811
  };
  int h[8] = {
    101527, 101531, 101533, 101537,
    101561, 101573, 101581, 101599
  };
    
  const char *ptr = (char*) &value;

  for (int i = 0; i < sizeof(unsigned int); ++i) {
    for (int j = 0; j < 8; ++j) {
      h[j] = (h[j] * p[j])^ptr[i];
    }
  }

  unsigned int ret = h[0];
  for (int i = 1; i < 8; ++i) {
    ret ^= h[i];
  }

  return h[0];
}
