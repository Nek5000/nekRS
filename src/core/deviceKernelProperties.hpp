#ifndef deviceKernelProperties_hpp_
#define deviceKernelProperties_hpp_

#include "nekrsSys.hpp"

class deviceKernelProperties {
 public:
  explicit deviceKernelProperties(occa::json properties)
      : properties_{properties} {}

  operator occa::json() const { return properties_; }
  operator occa::json&() { return properties_; }
  
  occa::json& define(const std::string& s) { 
    return properties_["defines"][s]; 
  }

  occa::json& include() { 
    return properties_["includes"].asArray(); 
  }

  occa::json& compiler_flags() { 
    return properties_["compiler_flags"].asArray(); 
  }

  occa::json& okl_include_paths() {
    return properties_["okl"]["include_paths"].asArray();
  }

 private:
  occa::json properties_;
};

#endif
