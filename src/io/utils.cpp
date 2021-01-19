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

void copyFile(const char *srcName, const char* destName)
{
  std::fstream src,dest;
  src.open (srcName);
  dest.open (destName);

  std::filebuf* inbuf  = src.rdbuf();
  std::filebuf* outbuf = dest.rdbuf();

  char c = inbuf->sbumpc();
  while (c != EOF)
  {
    outbuf->sputc (c);
    c = inbuf->sbumpc();
  }

  dest.close();
  src.close();
}
