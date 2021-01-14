#include "filesystem.hpp"

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 500 // Get nftw() declarations
#endif

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <ftw.h>
#include <sys/stat.h>
#include <limits.h>
#include <libgen.h>

static int rm_from_nftw(const char* pathname,
                        const struct stat* sbuf,
                        int type,
                        struct FTW* ftwb)
{
  switch (type) {
  case FTW_F:
  case FTW_DP:
    if (remove(pathname) != 0) {
      std::cerr << "ERROR: could not remove " << pathname << std::endl;
      return -1;
    }
    break;
  default:
    break;
  }
  return 0;
}

std::deque<std::string> fsys_path_split(const std::string& path)
{
  std::deque<std::string> result;
  char p[PATH_MAX];
  strcpy(p, path.c_str());

  do {
    auto dir = dirname(p);
    auto base = basename(p);
    result.push_front(base);
    strcpy(p, dir);
  } while (strcmp(p, ".") != 0 && strcmp(p, "/") != 0);

  if (strcmp(p, "/") == 0) {
    result.push_front("/");
  }

  return result;
}

int fsys_mkdir(const std::string& path, bool recurse)
{
  std::deque<std::string> dirs;
  if (recurse) {
    dirs = fsys_path_split(path);
  } else {
    dirs.push_back(path);
  }

  std::string next_dir;
  for (const auto& d : dirs) {
    next_dir = next_dir + "/" + d;
    struct stat sb;
    if (stat(next_dir.c_str(), &sb) == 0) {
      if (!S_ISDIR(sb.st_mode)) {
        std::cerr << "ERROR: cannot create directory '" << next_dir
                  << "': not a directory" << std::endl;
        return -1;
      }
    } else if (mkdir(next_dir.c_str(), S_IRWXU | S_IRWXG) != 0) {
      std::cerr << "ERROR: failed to create directory `" << next_dir << "'" << std::endl;
      return -1;
    }
  }

  return 0;
}

int fsys_rm(const std::string &path, bool followSymlinks)
{
  int flags = FTW_DEPTH;
  if (!followSymlinks) {
    flags |= FTW_PHYS;
  }

  struct stat s;
  if (stat(path.c_str(), &s) == 0) {
    if (nftw(path.c_str(), rm_from_nftw, 10, flags) != 0) {
      std::cerr << "ERROR: could not recursively remove " << path << std::endl;
      return -1;
    }
  }
  return 0;
}

