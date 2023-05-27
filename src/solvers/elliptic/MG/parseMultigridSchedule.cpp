#include "parseMultigridSchedule.hpp"
#include "nrs.hpp"

#include <set>
#include <limits>

std::pair<std::map<std::pair<int, bool>, int>, std::string>
parseMultigridSchedule(const std::string &schedule, setupAide& options, int defaultDegree)
{
  std::string errorString;
  std::map<std::pair<int, bool>, int> scheduleMap;
  auto entries = serializeString(schedule, ',');

  std::vector<int> levels;

  const auto INVALID = std::numeric_limits<int>::lowest();
  int prevOrder = std::numeric_limits<int>::max();
  int minOrder = std::numeric_limits<int>::max();

  // compute minimum order
  for (auto &&entry : entries) {
    auto tokens = serializeString(entry, '+');
    for (auto &&token : tokens) {
      if (token.find("p") != std::string::npos) {
        const auto order = std::stoi(serializeString(token, '=').at(1));
        minOrder = std::min(minOrder, order);
      }
    }
  }

  bool downLeg = true;

  for (auto &&entry : entries) {
    auto tokens = serializeString(entry, '+');

    int order = INVALID;
    int degree = INVALID;
    for (auto &&token : tokens) {
      if (token.find("p") != std::string::npos) {
        order = std::stoi(serializeString(token, '=').at(1));

        if (order > prevOrder) {
          downLeg = false;
        }

        if(downLeg) levels.push_back(order);
        
        prevOrder = order;
      }
      else if (token.find("degree") != std::string::npos) {
        degree = std::stoi(serializeString(token, '=').at(1));
      }
      else if (token.find("sweeps") != std::string::npos) {
        errorString += "setting number of sweeps is currently not supported!\n";
      }
      else {
        errorString += "unknown token '" + token + "' in schedule '" + schedule + "'!\n";
      }
    }

    if (order == INVALID) {
      errorString += "order not specified in " + entry + "\n";
    }

    if (degree == INVALID && order != minOrder){
      degree = defaultDegree;
    }

    scheduleMap[{order, downLeg}] = degree;
  }

  // up leg and down leg orders must be identical
  std::set<int> downLegOrders;
  std::set<int> upLegOrders;
  for (auto &&entry : scheduleMap) {
    const auto order = entry.first.first;
    if (entry.first.second) {
      if(order > minOrder){
        downLegOrders.insert(order);
      }
    }
    else {
      if(order > minOrder){
        upLegOrders.insert(order);
      }
    }
  }

  if (downLegOrders != upLegOrders) {
    errorString += "down leg and up leg orders must be identical!\n";
  }

  // all orders, except the coarse grid, must have a degree associated with them
  for (auto &&entry : scheduleMap) {
    if (entry.first.first == minOrder) {
      continue;
    }
    if (entry.second == INVALID) {
      errorString += "degree not specified for order " + std::to_string(entry.first.first) + "!\n";
    }
  }

  // all levels are positive
  for(auto&& level : levels){
    if(level < 1){
      errorString += "encountered level = " + std::to_string(level) + " < 1!\n";
    }
  }

  // top level order must match order of simulation
  int N;
  options.getArgs("POLYNOMIAL DEGREE", N);
  if(N != levels.front()){
    errorString += "finest degree " + std::to_string(levels.front())
      + " does not match polynomial degree " + std::to_string(N) + "!\n";
  }

#if 0
  // each successive level must be smaller
  for (unsigned i = 0U; i < levels.size(); ++i) {
    if (i > 0){
      if(levels[i] >= levels[i - 1]){
        errorString += "order[i] = " + std::to_string(levels[i]) + " >= order[i-1] = " + std::to_string(levels[i-1]) + "!\n";
        errorString += "\tEach level in the downward leg of the V-cycle must have order less than the previous level!\n";
      }
    }
  }
#endif

  return {scheduleMap, errorString};
}
