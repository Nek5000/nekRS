#include <cstdlib>
#include <fstream>
#include <cstdio>

#include <unistd.h>
#include <sys/stat.h>

#include "environment.hpp"

namespace os{
  std::string getWorkingDir(){
     char temp[FILENAME_MAX];
     return getcwd(temp,sizeof(temp)) ? std::string(temp) : std::string("");
  }

  std::string joinPath(std::string root,std::string relative){
    return root+separator_+relative;
  }

  int mkDir(std::string dirName){
    return mkdir(dirName.c_str(),S_IRWXU);
  }
}

namespace env{
  static int readConfig(){
    int rank; MPI_Comm_rank(comm_,&rank);
    if(rank==0){
      std::ifstream config(os::joinPath(shareDir(),"config"));
      std::string line;
      if(config.is_open()){
        while(std::getline(config,line))
          m_.insert(
            std::pair<std::string,std::string>(
              line.substr(0,line.find(delim_)),line.erase(0,line.find(delim_)+delim_.length())
            )
          );
      } else return 1;
    }
    return 0;
  }

  int init(MPI_Comm comm){
    MPI_Comm_dup(comm,&comm_);
    readConfig();
  }

  std::string installDir(){
    const char *dir =std::getenv("NEKRS_INSTALL_DIR");
    return dir ? std::string(dir) : installDir_;
  }

  std::string binDir(){
    const char *bin=std::getenv("NEKRS_BIN_DIR");
    return bin ? std::string(bin) : os::joinPath(installDir(),"bin");
  }

  std::string shareDir(){
    const char *share=std::getenv("NEKRS_SHARE_DIR");
    return share ? std::string(share) : os::joinPath(installDir(),"share");
  }

  std::string nekDir(){
    const char *nek=std::getenv("NEKRS_NEK5000_DIR");
    return nek ? std::string(nek) : os::joinPath(installDir(),"nek5000");
  }

  std::string libPDir(){
    const char *libP=std::getenv("NEKRS_LIBP_DIR");
    return libP ? std::string(libP) : os::joinPath(installDir(),"libparanumal");
  }

  std::string udfDir(){
    const char *udf=std::getenv("NEKRS_UDF_DIR");
    return udf ? std::string(udf) : os::joinPath(installDir(),"udf");
  }

  std::string nekInterfaceDir(){
    const char *nekInterface=std::getenv("NEKRS_NEKINTERFACE_DIR");
    return nekInterface ? std::string(nekInterface) : os::joinPath(installDir(),"nekInterface");
  }

  std::string occaDir(){
    const char *occa=std::getenv("OCCA_DIR");
    return occa ? std::string(occa) : os::joinPath(installDir(),"occa");
  }

  std::string cacheDir(){
    const char *cache=std::getenv("NEKRS_CACHE_DIR");
    if(cache==NULL){
      auto cacheDir_=os::joinPath(os::getWorkingDir(),".cache");
      setenv("NEKRS_CACHE_DIR",cacheDir_.c_str(),1);
      int rank; MPI_Comm_rank(comm_,&rank);
      if(rank==0) os::mkDir(cacheDir_);
      cache=cacheDir_.c_str();
    }
    return std::string(cache);
  }

  std::string getConfigValue(std::string key){
    const char *val=std::getenv(key.c_str());
    if(val==NULL){
      auto it=m_.find(key);
      if(it!=m_.end()) val=it->second.c_str();
    }
    return std::string(val);
  }

  std::string cxxCompiler(){ return getConfigValue("NEKRS_CXX"); }
  std::string cxxFlags(){ return getConfigValue("NEKRS_CXXFLAGS"); }
  std::string cCompiler(){ return getConfigValue("NEKRS_CC"); }
  std::string cFlags(){ return getConfigValue("NEKRS_CFLAGS"); }
  std::string fortranCompiler(){ return getConfigValue("NEKRS_FC"); }
  std::string fortranFlags(){ return getConfigValue("NEKRS_FFLAGS"); }
}
