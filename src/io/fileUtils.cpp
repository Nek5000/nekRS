#include <sys/stat.h>
#include <fstream> 
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <libgen.h>
#include <string>
#include <cstring>
#include <iostream>
#include <filesystem>

#include "device.hpp"
#include "platform.hpp"
#include "fileUtils.hpp"

void fileSync(const char * file)
{
  const std::string dir(fs::path(file).parent_path());

  int fd; 
  fd = open(file, O_RDONLY);
  fsync(fd);
  close(fd);

  fd = open(dir.c_str(), O_RDONLY);
  fsync(fd);
  close(fd);
}

bool _mkdir(const fs::path &file_path)
{
  size_t pos = 0;
  auto ret_val = true;

  std::string dir_path(file_path);
  if (!fs::is_directory((file_path)))
    dir_path = file_path.parent_path();

  while (ret_val && pos != std::string::npos) {
    pos = dir_path.find('/', pos + 1);
    const auto dir = fs::path(dir_path.substr(0, pos));
    if (!fs::exists(dir)) {
      ret_val = fs::create_directory(dir);
    }
  }

  return ret_val;
}

void fileBcast(const fs::path &srcPathIn,
               const fs::path &dstPath,
               MPI_Comm comm,
               int verbose)
{
  int rank;
  MPI_Comm_rank(comm, &rank);

  const auto path0 = fs::current_path();

  int localRank;
  const int localRankRoot = 0;
  MPI_Comm commLocal;
  {
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &commLocal);
    MPI_Comm_rank(commLocal, &localRank);
  }

  int nodeRank = -1;
  const int nodeRankRoot = 0;
  MPI_Comm commNode;
  {
    int nodeColor = MPI_UNDEFINED;
    if (localRank == localRankRoot)
      nodeColor = 1;

    MPI_Comm_split(comm, nodeColor, rank, &commNode);
    if (commNode != MPI_COMM_NULL)
      MPI_Comm_rank(commNode, &nodeRank);
  }


  std::vector<std::string> fileList;
  if (nodeRank == nodeRankRoot) {

    nrsCheck(!fs::exists(srcPathIn), MPI_COMM_SELF, EXIT_FAILURE, 
             "Cannot find %s!\n", std::string(srcPathIn).c_str());

    const auto srcPathCanonical = fs::canonical(srcPathIn);
    fs::current_path(srcPathCanonical.parent_path()); 
    const auto srcPath = fs::relative(srcPathCanonical, fs::current_path());

    if (!fs::is_directory((srcPath))) {
      fileList.push_back(srcPath);
    } else {
      for (const auto &entry : fs::recursive_directory_iterator(srcPath)) {
        if (entry.is_regular_file() || entry.is_symlink()) {
          fileList.push_back(entry.path());
        }
      }
    }

  }
  int maxFileSize = 0;
  int nFiles = (nodeRank == nodeRankRoot) ? fileList.size() : 0;
  MPI_Bcast(&nFiles, 1, MPI_INT, nodeRankRoot, comm);

  for (int i = 0; i < nFiles; i++) {
    int bufSize = (nodeRank == nodeRankRoot) ? fileList.at(i).size() : 0;
    MPI_Bcast(&bufSize, 1, MPI_INT, nodeRankRoot, comm);

    auto buf = (char *)std::malloc(bufSize * sizeof(char));
    if (nodeRank == nodeRankRoot) {
      std::strncpy(buf, fileList.at(i).c_str(), bufSize);
      auto fileSize = fs::file_size(fileList.at(i));
      nrsCheck(fileSize > std::numeric_limits<int>::max(),
               MPI_COMM_SELF, EXIT_FAILURE,
               "%s\n", "File size too large!");
      maxFileSize = std::max(maxFileSize, static_cast<int>(fileSize));
    }
    MPI_Bcast(buf, bufSize, MPI_CHAR, nodeRankRoot, comm);
    if (nodeRank != nodeRankRoot)
      fileList.push_back(std::string(buf, 0, bufSize));
    free(buf);
  }

  MPI_Bcast(&maxFileSize, 1, MPI_INT, nodeRankRoot, comm);

  auto fileBuf = (char *) std::malloc(maxFileSize * sizeof(char)); 
  nrsCheck(fileBuf == nullptr, MPI_COMM_SELF, 
           EXIT_FAILURE, "%s\n", "allocating file buffer failed!");

  // sweep through list and transfer to nodes 
  for (const auto &file : fileList) {
    int bufSize = 0;
    const std::string filePath = dstPath / fs::path(file);

    if (commNode != MPI_COMM_NULL) {
      if (nodeRank == nodeRankRoot)
        bufSize = fs::file_size(file);
      MPI_Bcast(&bufSize, 1, MPI_INT, nodeRankRoot, commNode);

      if (nodeRank == nodeRankRoot) {
        std::ifstream input(file, std::ifstream::binary);
        input.read(fileBuf, bufSize);
        input.close();
      }
      MPI_Bcast(fileBuf, bufSize, MPI_CHAR, nodeRankRoot, commNode);

      if (nodeRank == nodeRankRoot && verbose)
        std::cout << __func__ << ": " << file << " -> " << std::flush;
    }

    // write file collectively to node-local storage;
    if (localRank == localRankRoot) {
      _mkdir(filePath);
      fs::remove(filePath); 
    }

#if 0
    MPI_File fh;
    int retVal = MPI_File_open(commLocal, filePath.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    if (retVal) {
      char err_string[MPI_MAX_ERROR_STRING];
      int err_string_len;
      MPI_Error_string(retVal, err_string, &err_string_len);
      nrsAbort(MPI_COMM_SELF, EXIT_FAILURE,
               "MPI_File_open error: %s\n", err_string);
    }

    if (localRank == localRankRoot) {
      MPI_Status status;
      retVal = MPI_File_write_at(fh, 0, fileBuf, bufSize, MPI_BYTE, &status);
      nrsCheck(retVal, MPI_COMM_SELF, EXIT_FAILURE,
               "MPI_File_write_at with retVal=%d\n!", retVal);
    }

    MPI_File_sync(fh);
    MPI_Barrier(commLocal);
    MPI_File_sync(fh);
    retVal = MPI_File_close(&fh);
    nrsCheck(retVal, MPI_COMM_SELF, EXIT_FAILURE,
               "MPI_File_close with retVal=%d\n!", retVal);
#else
    if (localRank == localRankRoot) {
      std::ofstream fh (filePath, std::ios::out | std::ios::binary | std::ios::trunc);
      fh.write(fileBuf, bufSize);
      fh.close();
      fileSync(filePath.c_str());
    }
    MPI_Barrier(commLocal);
#endif

    if (localRank == localRankRoot) {
      fs::permissions(filePath, fs::perms::owner_all);
    }

    if (commNode != MPI_COMM_NULL) {
      if (nodeRank == nodeRankRoot && verbose) {
        std::cout << fs::canonical(filePath) << " (" << bufSize << " bytes)" << std::endl;
      }
    }
  }

  free(fileBuf);
  MPI_Comm_free(&commLocal);
  if(commNode != MPI_COMM_NULL) MPI_Comm_free(&commNode);
  fs::current_path(path0);
  MPI_Barrier(comm);
}

bool isFileNewer(const char *file1, const char *file2)
{
  struct stat s1, s2;
  lstat(file1, &s1);
  if (lstat(file2, &s2) != 0) return true; 
  if (s1.st_mtime > s2.st_mtime) 
    return true;
  else
    return false;	  
}

void copyFile(const char *srcFile, const char *dstFile)
{
  std::ifstream src(srcFile, std::ios::binary);
  std::ofstream dst(dstFile, std::ios::trunc | std::ios::binary);
  dst << src.rdbuf();
  src.close();
  dst.close();
  fileSync(dstFile);
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
