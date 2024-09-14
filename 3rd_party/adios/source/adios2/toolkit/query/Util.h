#ifndef ADIOS2_QUERY_UTIL_H
#define ADIOS2_QUERY_UTIL_H

#include <cctype>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <stdexcept>
#include <string>
#include <vector>

namespace adios2
{
namespace query
{
/*
static size_t ToUIntValue(const adios2::Params &params, const std::string &key,
                        size_t defaultValue)
{
  auto it = params.find(key);
  if (it != params.end())
  {
      try
      {
          auto value = (size_t)(std::stoul(it->second));
          return value;
      }
      catch (std::exception &e)
      {
          std::cout << e.what() << std::endl;
          return defaultValue;
      }
  }
  return defaultValue;
}
*/

// The next three are implemented in Worker.cpp
bool EndsWith(const std::string &hostStr, const std::string &fileTag);
bool IsFileNameXML(const std::string &filename);
bool IsFileNameJSON(const std::string &filename);

}; // namespace query
}; // name space adios2

#endif // QUERY_WORKER_H
