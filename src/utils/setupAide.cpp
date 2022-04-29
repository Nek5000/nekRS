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

#include <tuple>
#include "setupAide.hpp"

std::string setupAide::getArgs(std::string key) const
{
  if(keyWordToDataMap.count(key) == 0){
    return "";
  }
  return keyWordToDataMap.at(key);
}

void setupAide::removeArgs(std::string key)
{
  auto iter = keyWordToDataMap.find(key);
  if(iter != keyWordToDataMap.end()){
    keyWordToDataMap.erase(iter);
  }
}

void setupAide::setArgs(std::string key, std::string value)
{
  keyWordToDataMap[key] = value;
}

int setupAide::getArgs(std::string key, std::vector < std::string >& m, std::string delimeter) const
{
  std::string args, current;
  std::vector < std::string > argv;
  int argc, size;

  args = getArgs(key);

  size = args.length();

  current = "";

  for(int i = 0; i < size; i++) { // TW
    while( i < size && delimeter.find(args[i]) == std::string::npos )
      current += args[i++];

    if(current.length())
      argv.push_back(current);

    current = "";
  }

  argc = argv.size();

  if(!argc)
    return 0;

  m.resize(argc);

  for(int i = 0; i < argc; i++) // TW
    m[i] = argv[i];

  return 1;
}

int setupAide::compareArgs(std::string key, std::string token) const
{
  std::string foundToken;
  if(getArgs(key,foundToken)) {
    if(foundToken == token)
      return 1;
    if(foundToken.find(token) != std::string::npos)
      return 2;
  }

  return 0;
}

std::ostream & operator << (std::ostream &os, const setupAide &aide){
  int maxLength = 0;
  for(auto&& keyAndValuePair : aide.keyWordToDataMap)
  {
    const std::string key = keyAndValuePair.first;
    int L = key.length();
    if(L > maxLength)
      maxLength = L;
  }

  std::string key, value;

  for(auto&& keyAndValuePair : aide.keyWordToDataMap)
  {
    std::tie(key, value) = keyAndValuePair;
    os << "key: " << key << ",";

    for(int j = key.length(); j < maxLength; ++j)
      os << " ";

    os << "value: " << value << std::endl;
  }

  return os;
  }
