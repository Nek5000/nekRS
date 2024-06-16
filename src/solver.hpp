#if !defined(nekrs_solve_hpp_)
#define nekrs_solve_hpp_

#include "nekrsSys.hpp"
#include "mesh3D.h"

class solver_t {
  public:
    solver_t() {}; 
    virtual ~solver_t() {};

    virtual std::string id() const { return ""; };

    virtual void printMinMax() {};
    virtual void printRunStat(int step) {};
};

#endif
