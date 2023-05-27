#ifndef parse_multigrid_schedule_hpp_
#define parse_multigrid_schedule_hpp_

#include <string>
#include <vector>
#include <map>

class setupAide;

std::pair<std::map<std::pair<int, bool>, int>, std::string>
parseMultigridSchedule(const std::string &schedule, setupAide& options, int defaultOrder);

#endif