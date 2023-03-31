#include "parseMultigridSchedule.hpp"
#include "determineMGLevels.hpp"
#include "platform.hpp"

std::vector<int> determineMGLevels(std::string section)
{
  const std::string optionsPrefix = [section]() {
    std::string prefix = section + std::string(" ");
    if (section.find("temperature") != std::string::npos) {
      prefix = std::string("scalar00 ");
    }
    std::transform(prefix.begin(), prefix.end(), prefix.begin(), [](unsigned char c) {
      return std::toupper(c);
    });
    return prefix;
  }();

  std::vector<int> levels;
  int N;
  platform->options.getArgs("POLYNOMIAL DEGREE", N);

  std::string p_mgschedule = platform->options.getArgs(optionsPrefix + "MULTIGRID SCHEDULE");
  if(!p_mgschedule.empty()){

    // note: default order is not required here.
    // We just need the levels, not the degree.
    auto [scheduleMap, errorString] = parseMultigridSchedule(p_mgschedule, platform->options, 3);

    nrsCheck(errorString.size(), platform->comm.mpiComm, EXIT_FAILURE,
             "%s\n", errorString.c_str());

    for(auto && [cyclePosition, smootherOrder] : scheduleMap){
      auto [order, isDownLeg] = cyclePosition;
      if(isDownLeg){
        levels.push_back(order);
      }
    }

    std::sort(levels.rbegin(), levels.rend());

    if (levels.back() > 1) {
      if (platform->options.compareArgs(optionsPrefix + "MULTIGRID COARSE SOLVE", "TRUE")) {
        // if the coarse level has p > 1 and requires solving the coarsest level,
        // rather than just smoothing, SEMFEM must be used for the discretization
        const auto usesSEMFEM =
            platform->options.compareArgs(optionsPrefix + "MULTIGRID SEMFEM", "TRUE");

        nrsCheck(!usesSEMFEM, platform->comm.mpiComm, EXIT_FAILURE, 
                 "FEM coarse discretization only supports p=1 for the coarsest level!\n", "");
      }
    }

    return levels;
  }

  if (platform->options.compareArgs(optionsPrefix + "MULTIGRID SMOOTHER", "ASM") ||
           platform->options.compareArgs(optionsPrefix + "MULTIGRID SMOOTHER", "RAS")) {
    std::map<int, std::vector<int>> mg_level_lookup = {
        {1, {1}},
        {2, {2, 1}},
        {3, {3, 1}},
        {4, {4, 2, 1}},
        {5, {5, 3, 1}},
        {6, {6, 3, 1}},
        {7, {7, 3, 1}},
        {8, {8, 5, 1}},
        {9, {9, 5, 1}},
        {10, {10, 6, 1}},
        {11, {11, 6, 1}},
        {12, {12, 7, 1}},
        {13, {13, 7, 1}},
        {14, {14, 8, 1}},
        {15, {15, 9, 1}},
    };

    return mg_level_lookup.at(N);
  }

  std::map<int, std::vector<int>> mg_level_lookup = {
      {1, {1}},
      {2, {2, 1}},
      {3, {3, 1}},
      {4, {4, 2, 1}},
      {5, {5, 3, 1}},
      {6, {6, 4, 2, 1}},
      {7, {7, 5, 3, 1}},
      {8, {8, 6, 4, 1}},
      {9, {9, 7, 5, 1}},
      {10, {10, 8, 5, 1}},
      {11, {11, 9, 5, 1}},
      {12, {12, 10, 5, 1}},
      {13, {13, 11, 5, 1}},
      {14, {14, 12, 5, 1}},
      {15, {15, 13, 5, 1}},
  };

  return mg_level_lookup.at(N);
}
