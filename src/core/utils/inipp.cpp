/*
   MIT License

   Copyright (c) 2017-2018 Matthias C. M. Troffaes

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#include <cstring>
#include <string>
#include <iostream>
#include <list>
#include <map>
#include <algorithm>
#include <functional>
#include <cctype>
#include <sstream>
#include <vector>

#include <inipp.hpp>

namespace
{

static bool enforceLowerCase = false;

static std::vector<std::string> nothing = {};
static std::vector<std::string> generalKeys = {
  {"dt"},
  {"endTime"},
  {"numSteps"},
  {"polynomialOrder"},
  {"dealiasing"},
  {"cubaturePolynomialOrder"},
  {"startFrom"},
  {"stopAt"},
  {"elapsedtime"},
  {"timestepper"},
  {"subCyclingSteps"},
  {"subCycling"},
  {"writeControl"},
  {"writeInterval"},
  {"constFlowRate"},
  {"verbose"},
  {"variableDT"},

  {"oudf"},
  {"udf"},
  {"usr"},

};

static std::vector<std::string> problemTypeKeys = {
  {"stressFormulation"},
  {"equation"},
};

// common keys
static std::vector<std::string> commonKeys = {
  {"solver"},
  {"residualTol"},
  {"initialGuess"},
  {"preconditioner"},
  {"pMultigridCoarsening"},
  {"smootherType"},
  {"coarseSolver"},
  {"boundaryTypeMap"},
  {"maxIterations"},
  {"regularization"},

  // deprecated filter params
  {"filtering"},
  {"filterWeight"},
  {"filterModes"},
  {"filterCutoffRatio"},

  // deprecated no-op extrapolation param
  {"extrapolation"},


  // deprecated projection params
  {"residualProj"},
  {"residualProjection"},
  {"residualProjectionVectors"},
  {"residualProjectionStart"},
};

static std::vector<std::string> meshKeys = {
  {"partitioner"},
  {"file"},
};

static std::vector<std::string> velocityKeys = {
  {"density"},
  {"viscosity"},
};

static std::vector<std::string> temperatureKeys = {
  {"rhoCp"},
  {"conductivity"},
};

static std::vector<std::string> scalarKeys = {
  {"rho"},
  {"diffusivity"},
};

static std::vector<std::string> boomeramgKeys = {
  {"coarsenType"},
  {"interpolationType"},
  {"smootherType"},
  {"iterations"},
  {"strongThreshold"},
  {"nonGalerkinTol"},
  {"aggressiveCoarseningLevels"},
};

static std::vector<std::string> amgxKeys = {
  {"configFile"},
};
static std::vector<std::string> occaKeys = {
  {"backend"},
  {"deviceNumber"},
};

static std::vector<std::string> pressureKeys = {};

static std::vector<std::string> deprecatedKeys = {
  // deprecated filter params
  {"filtering"},
  {"filterWeight"},
  {"filterModes"},
  {"filterCutoffRatio"},

  // deprecated no-op extrapolation param
  {"extrapolation"},

  // deprecated projection params
  {"residualProj"},
  {"residualProjection"},
  {"residualProjectionVectors"},
  {"residualProjectionStart"},
};

void convertToLowerCase(std::vector<std::string>& stringVec)
{
  for(auto && s : stringVec){
    std::transform(s.begin(), s.end(), s.begin(),
      [](unsigned char c){ return std::tolower(c); });
  }
}

void makeStringsLowerCase()
{
  convertToLowerCase(generalKeys);
  convertToLowerCase(problemTypeKeys);
  convertToLowerCase(commonKeys);
  convertToLowerCase(meshKeys);
  convertToLowerCase(temperatureKeys);
  convertToLowerCase(scalarKeys);
  convertToLowerCase(deprecatedKeys);
  convertToLowerCase(amgxKeys);
  convertToLowerCase(boomeramgKeys);
  convertToLowerCase(pressureKeys);
  convertToLowerCase(occaKeys);
}

const std::vector<std::string>& getValidKeys(const std::string& section)
{
  if(!enforceLowerCase)
  {
    makeStringsLowerCase();
    enforceLowerCase = true;
  }

  if(section == "general")
    return generalKeys;
  if(section == "problemtype")
    return problemTypeKeys;
  if(section == "mesh")
    return meshKeys;
  if(section == "temperature")
    return temperatureKeys;
  if(section == "pressure")
    return pressureKeys;
  if(section.find("scalar") != std::string::npos)
    return scalarKeys;
  if(section == "amgx")
    return amgxKeys;
  if(section == "boomeramg")
    return boomeramgKeys;
  if(section == "occa")
    return occaKeys;
  if(section == "velocity")
    return velocityKeys;
  else
    return nothing;
}
}


namespace inipp
{
namespace detail
{

inline void ltrim(std::string & s)
{
  s.erase(s.begin(),
          std::find_if(s.begin(), s.end(),
                       [](int ch) {
          return !std::isspace(ch);
        }));
}

inline void rtrim(std::string & s)
{
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](int ch) {
          return !std::isspace(ch);
        }).base(),
          s.end());
}

inline bool replace(std::string & str,
                    const std::string & from,
                    const std::string & to)
{
  auto changed = false;
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
    changed = true;
  }
  return changed;
}
} // namespace detail

bool Ini::extract(const std::string & key,
             const std::string & value,
             std::string & dst)
{
  if (sections[key].count(value)) {
    dst = sections[key][value];
    return true;
  } else {
    return false;
  }
}

void Ini::generate(std::ostringstream & os) const
{
  for (auto const & sec : sections) {
    os << char_section_start << sec.first << char_section_end << std::endl;
    for (auto const & val : sec.second)
      os << val.first << char_assign << val.second << std::endl;
  }
}

void Ini::parse(std::stringstream & is, bool lowerValue)
{
  std::string line;
  std::string section;
  while (!is.eof()) {
    std::getline(is, line);
    auto it = std::find_if(line.rbegin(), line.rend(),
                           [](int ch) { return ch == char_comment; });
    if (it != line.rend()) line.erase((++it).base(), line.end());
    detail::ltrim(line);
    detail::rtrim(line);

    const auto length = line.length();
    if (length > 0) {
      const auto pos = line.find_first_of(char_assign);
      const auto & front = line.front();
      if (front == char_comment) {
        continue;
      }else if (front == char_section_start) {
        if (line.back() == char_section_end) {
          section = line.substr(1, length - 2);
          transform(section.begin(), section.end(), section.begin(),
                    std::ptr_fun<int, int>(std::tolower));
        } else {
          errors.push_back(line);
        }
      }else if (pos != 0 && pos != std::string::npos) {
        std::string variable(line.substr(0, pos));
        std::string value(line.substr(pos + 1, length));
        transform(variable.begin(), variable.end(), variable.begin(),
                  std::ptr_fun<int, int>(std::tolower));
        detail::rtrim(variable);
        detail::ltrim(value);

        bool inquotes = lowerValue && value.front() == '"' &&
                        lowerValue && value.back() == '"';
        if (lowerValue && !inquotes)
          transform(value.begin(), value.end(), value.begin(),
                    std::ptr_fun<int, int>(std::tolower));
        value.erase(std::remove(value.begin(), value.end(), '"'), value.end());

        auto & sec = sections[section];
        if (sec.find(variable) == sec.end())
          sec.insert(std::make_pair(variable, value));
        else
          errors.push_back(line);
      }else {
        errors.push_back(line);
      }
    }
  }
}

int Ini::validateKeys() const
{
  int err = 0;
  for (auto const & sec : sections) {
    if(sec.first.find("caseparams") != std::string::npos) continue;
    const auto& validKeys = getValidKeys(sec.first);
    for (auto const & val : sec.second) {
      if (std::find(validKeys.begin(), validKeys.end(), val.first) == validKeys.end()) {
        if (std::find(commonKeys.begin(), commonKeys.end(), val.first) == commonKeys.end()) {
          std::cout << "par-file: " << sec.first << "." << val.first << " unknown!\n";
          err++;
        }
      }
    }
  }
  return err;
}

void Ini::printDeprecation() const
{
  for (auto const & sec : sections) {
    for (auto const & val : sec.second) {
      if (std::find(deprecatedKeys.begin(), deprecatedKeys.end(), val.first) != deprecatedKeys.end()) {
          std::cout << "par-file: " << sec.first << "." << val.first << " deprecated!\n";
      }
    }
  }
}

void Ini::interpolate()
{
  int global_iteration = 0;
  auto changed = false;
  // replace each "${variable}" by "${section:variable}"
  for (auto & sec : sections)
    replace_symbols(local_symbols(sec.first, sec.second), sec.second);
  // replace each "${section:variable}" by its value
  do {
    changed = false;
    const auto syms = global_symbols();
    for (auto & sec : sections)
      changed |= replace_symbols(syms, sec.second);
  } while (changed && (max_interpolation_depth > global_iteration++));
}

void Ini::default_section(const Ini::Section & sec)
{
  for (auto & sec2 : sections)
    for (const auto & val : sec)
      sec2.second.insert(val);
}

void Ini::clear()
{
  sections.clear();
  errors.clear();
}

std::string Ini::local_symbol(const std::string & name) const
{
  return char_interpol + (char_interpol_start + name + char_interpol_end);
}

std::string Ini::global_symbol(const std::string & sec_name, const std::string & name) const
{
  return local_symbol(sec_name + char_interpol_sep + name);
}

Ini::Symbols Ini::local_symbols(const std::string & sec_name, const Ini::Section & sec) const
{
  Ini::Symbols result;
  for (const auto & val : sec)
    result.push_back(std::make_pair(local_symbol(val.first), global_symbol(sec_name, val.first)));
  return result;
}

Ini::Symbols Ini::global_symbols() const
{
  Ini::Symbols result;
  for (const auto & sec : sections)
    for (const auto & val : sec.second)
      result.push_back(
        std::make_pair(global_symbol(sec.first, val.first), val.second));
  return result;
}

bool Ini::replace_symbols(const Ini::Symbols & syms, Ini::Section & sec) const
{
  auto changed = false;
  for (auto & sym : syms)
    for (auto & val : sec)
      changed |= detail::replace(val.second, sym.first, sym.second);
  return changed;
}

} // namespace inipp
