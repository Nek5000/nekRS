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

#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols){

  char buf[BUFSIZ];
  rewind(fp); // rewind to beginning

  int flag = 0;
  int status;
  char *stat;

  //search for label in file
  while(fgets(buf, BUFSIZ, fp)){
    if (strstr(buf, label)) {flag =1; break;};
  }
  if (flag==0) {
    printf("ERROR: Unable to find label: '%s' in node file.\n", label);
    exit(-1);
  }

  //if found read in number of rows and columns
  status = fscanf(fp, "%d %d",  Nrows, Ncols);
  stat = fgets(buf, BUFSIZ, fp); //read to end of line

  *A = (dfloat*) calloc((*Nrows)*(*Ncols), sizeof(dfloat)); //allocate space
  for(int n=0;n<(*Nrows)*(*Ncols);++n) //read matrix data
    status = fscanf(fp, dfloatFormat, (*A)+n);
}

void readIntArray(FILE *fp, const char *label, int **A, int *Nrows, int* Ncols){

  char buf[BUFSIZ];
  rewind(fp); // rewind to beginning

  int flag = 0;
  int status;
  char *stat;

  //search for label in file
  while(fgets(buf, BUFSIZ, fp)){
    if (strstr(buf, label)) {flag =1; break;};
  }
  if (flag==0) {
    printf("ERROR: Unable to find label: '%s' in node file.\n", label);
    exit(-1);
  }

  //if found read in number of rows and columns
  status = fscanf(fp, "%d %d",  Nrows, Ncols);
  stat = fgets(buf, BUFSIZ, fp); //read to end of line

  *A = (int*) calloc((*Nrows)*(*Ncols), sizeof(int)); //allocate space
  for(int n=0;n<(*Nrows)*(*Ncols);++n) //read matrix data
    status = fscanf(fp, "%d", (*A)+n);
}
