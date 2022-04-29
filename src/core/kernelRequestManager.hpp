#ifndef kernelRequestManager_hpp_
#define kernelRequestManager_hpp_
#include <occa.hpp>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <iostream>

class platform_t;

class kernelRequestManager_t
{

  struct kernelRequest_t
  {
    inline bool operator==(const kernelRequest_t& other) const
    {
      return requestName == other.requestName;
    }
    inline bool operator<(const kernelRequest_t& other) const
    {
      return requestName < other.requestName;
    }
    inline bool operator> (const kernelRequest_t& other) const { return *this < other; }
    inline bool operator<=(const kernelRequest_t& other) const { return !(*this > other); }
    inline bool operator>=(const kernelRequest_t& other) const { return !(*this < other); }
    inline bool operator!=(const kernelRequest_t& other) const { return !(*this == other); }

    kernelRequest_t(const std::string& m_requestName,
                    const std::string& m_fileName,
                    const occa::properties& m_props,
                    std::string m_suffix = std::string())
    :
    requestName(m_requestName),
    fileName(m_fileName),
    suffix(m_suffix),
    props(m_props)
    {}
    std::string requestName;
    std::string fileName;
    std::string suffix;
    occa::properties props;

    std::string to_string() const {
      std::ostringstream ss;
      ss << "requestName : " << requestName << "\n";
      ss << "fileName : " << fileName << "\n";
      ss << "suffix : " << suffix << "\n";
      ss << "props : " << props << "\n";;
      return ss.str();
    }
  };
public:
  kernelRequestManager_t(const platform_t& m_platform);
  void add(const std::string& m_requestName,
                  const std::string& m_fileName,
                  const occa::properties& m_props,
                  std::string m_suffix = std::string(),
                  bool assertUnique = false);
  void add(const std::string& requestName, occa::kernel kernel){
    requestToKernelMap[requestName] = kernel;
  }
  
  void compile();

  occa::kernel
  get(const std::string& request, bool checkValid = true) const;

  bool
  processed() const { return kernelsProcessed; }

private:
  const platform_t& platformRef;
  bool kernelsProcessed;
  std::set<kernelRequest_t> kernels;
  std::map<std::string, occa::kernel> requestToKernelMap;
  std::map<std::string, std::set<kernelRequest_t>> fileNameToRequestMap;

  void add(kernelRequest_t request, bool assertUnique = true);

};
#endif /** kernelRequestManager_hpp_ **/