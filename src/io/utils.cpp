#include <sys/stat.h>
#include <assert.h>
#include <fstream> 
#include <cstdio>

bool isFileNewer(const char *file1, const char* file2)
{
  struct stat s1, s2;
  if (lstat(file1, &s1) != 0) assert(1);
  if (lstat(file2, &s2) != 0) return true; 
  if (s1.st_mtime > s2.st_mtime) 
    return true;
  else
    return false;	  
}

void copyFile(const char *srcFile, const char* dstFile)
{
  std::ifstream src (srcFile, std::fstream::binary);
  std::ofstream dst (dstFile, std::fstream::trunc|std::fstream::binary);
  dst<<src.rdbuf();
  src.close();
  dst.close();
}

bool fileExists(const char *file)
{
  return realpath(file, NULL);
}

bool isFileEmpty(const char *file)
{
  std::ifstream f(file);
  const bool isEmpty = f.peek() == std::ifstream::traits_type::eof();
  f.close();
  return isEmpty;
}
