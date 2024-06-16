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
public:
  kernelRequestManager_t(const platform_t& m_platform);
  void add(const std::string& m_requestName,
                  const std::string& m_fileName,
                  const occa::properties& m_props,
                  std::string m_suffix = std::string(),
                  bool checkUnique = false);

  void add(const std::string& requestName, occa::kernel kernel);
  
  void compile();

  occa::kernel load(const std::string& request, const std::string& kernelName = std::string());

  bool processed() const { return kernelsProcessed; }

private:
  struct kernelRequest_t
  {
    std::string requestName; // assumed to be unique
    std::string fileName;
    std::string suffix;
    occa::properties props;
    occa::kernel kernel;

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

    std::string to_string() const {
      std::ostringstream ss;
      ss << "requestName : " << requestName << "\n";
      ss << "fileName : " << fileName << "\n";
      ss << "suffix : " << suffix << "\n";
      ss << "props : " << props << "\n";
      return ss.str();
    }
  };

  const platform_t& platformRef;
  bool kernelsProcessed;
  std::set<kernelRequest_t> requests;

  // request name to request
  std::map<std::string, kernelRequest_t> requestMap;

  // request and kernel name to kernel
  std::map<std::tuple<kernelRequest_t, std::string>, occa::kernel> kernelMap;

  void add(kernelRequest_t request, bool assertUnique = true);

};

#endif
