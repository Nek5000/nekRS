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

//#include "headers2d.hpp"

#include "setupAide.hpp"

using std::stringstream;
using std::fstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;

setupAide::setupAide(){}

setupAide::setupAide(string setupFile){
  read(setupFile);
}

setupAide::setupAide(const setupAide& sa){
  *this = sa;
}

setupAide& setupAide::operator = (const setupAide& sa){
  int size = sa.data.size();

  data.resize(size);
  keyword.resize(size);

  for(int i=0; i<size; i++){ // TW
    data[i]    = sa.data[i];
    keyword[i] = sa.keyword[i];
  }

  return *this;

}

string setupAide::readFile(string filename){
  struct stat statbuf;

  FILE *fh = fopen(filename.c_str(), "r");
  if (fh == 0){
    printf("Failed to open: %s\n", filename.c_str());
    exit(1);
  }

  stat(filename.c_str(), &statbuf);
  char *source = (char *) malloc(statbuf.st_size + 1);
  size_t status = fread(source, statbuf.st_size, 1, fh);
  source[statbuf.st_size] = '\0';

  string ret = source;

  return ret;
}

void setupAide::read(string setupFile){
  vector<string> data2;
  vector<string> keyword2;

  string args = readFile(setupFile);

  int size = args.length();
  string current = "";
  stringstream ss;
  char c;

  for(int i=0; i<size; i++){
    c = args[i];

    // Batch strings together
    if(c == '\'' || c == '"'){
      current += c;
      i++;

      while(i < size && args[i] != c)
        current += args[i++];

      if(i >= size)
        break;

      if( i < (size-1) )
        current += args[i];
    }

    // Batch comments
    else if(c == '/' && i < size && args[i+1] == '*'){
      i += 2;

      while( args[i] != '*' || (i < size && args[i+1] != '/') )
        i++;

      if(i >= size)
        break;

      i++;
    }

    // Removing # comments
    else if(c == '#'){
      i++;

      while(i < size && args[i] != '\n')
        i++;
    }

    // Change \[\] to []
    else if(c == '\\' && i < size && (args[i+1] == '[' || args[i+1] == ']')){
      current += args[i+1];
      i += 2;
    }

    // Split keywords []
    else if(c == '['){
      data2.push_back(current);
      current = "";
      i++;

      while(i < size && args[i] != ']')
        current += args[i++];

      keyword2.push_back(current);
      current = "";
    }

    // Else add the character
    else 
      if(!isspace(c)) // new check to remove whitespace
        current += c;

    if(i >= (size-1) && current.length())
      data2.push_back(current);

  }

  int argc = (data2.size() - 1);

  data.resize(argc);
  keyword.resize(argc);

  for(int i=0; i<argc; i++){ // TW
    data[i]    = data2[i+1];
    keyword[i] = keyword2[i];
  }
}

string setupAide::getArgs(string key){

  for(int i=0; i<keyword.size(); i++) // TW
    if(!( keyword[i].compare(key) ))
      return data[i];

  //printf("Warning: Failed to find [%s].\n", key.c_str());
  return "";
}

void setupAide::setArgs(string key, string value){
  for(int i=0; i<keyword.size(); i++) // TW
    if(!( keyword[i].compare(key) )) {
      data[i] = value;
      return;
    }

  //add the key value pair
  keyword.push_back(key);
  data.push_back(value);
  return;
}

int setupAide::getArgs(string key, vector<string>& m, string delimeter){
  string args, current;
  vector<string> argv;
  int argc, size;

  args = getArgs(key);

  size = args.length();

  current = "";

  for(int i=0; i<size; i++){ // TW
    while( i < size && delimeter.find(args[i]) == string::npos )
      current += args[i++];

    if(current.length())
      argv.push_back(current);

    current = "";
  }

  argc = argv.size();

  if(!argc){
    //printf("Warning: Failed to find [%s].\n", key.c_str());
    return 0;
  }

  m.resize(argc);

  for(int i=0; i<argc; i++)// TW
    m[i] = argv[i];

  return 1;
}


int setupAide::compareArgs(string key, string token){

  string foundToken;
  if(getArgs(key,foundToken)){
    if(foundToken==token)
      return 1;
    if(foundToken.find(token) != string::npos)
      return 2;
  }

  return 0;
  
}

std::ostream & operator << (std::ostream &os, const setupAide &aide){
  int maxLength = 0;
  for(int i=0; i<aide.keyword.size(); i++){
    int L = aide.keyword[i].length();
    if(L>maxLength)
      maxLength = L;
  }
  for(int i=0; i<aide.keyword.size(); i++){
    os << "key: " << aide.keyword[i] << ",";
    
    for(int j=aide.keyword[i].length();j<maxLength;++j)
      os << " ";
    
    os << "value: " << aide.data[i] << std::endl;
  }
  
  return os;
}
