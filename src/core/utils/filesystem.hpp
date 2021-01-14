#ifndef NEKRS_FILESYSTEM_H
#define NEKRS_FILESYSTEM_H

#include <string>
#include <deque>

// Split a given path into components
std::deque<std::string> fsys_path_split(const std::string& path);

// Make a directory at a given path
int fsys_mkdir(const std::string& path, bool recurse = true);

// Removes a file or directory.  
// * Always recursive and does not raise errors if a given file/dir 
//   does not exist (like rm -rf)
// * By default, will not follow symlinked dirs and will just 
//   remove the the symlinks themselves.  Set `followSymlinks` to 
//   `true` to follow symlinked dirs to the destination path and remove 
//   the contents at the destination.  
int fsys_rm(const std::string &path, bool followSymlinks = false);

#endif 
