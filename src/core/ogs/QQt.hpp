#ifndef QQT_HPP
#define QQT_HPP

#include "nekrsSys.hpp"
#include "ogs.hpp"


class QQt 
{

public:

struct timingConfig
{
  dlong stride = 0;
  int nVec = 1;
  std::string type = ogsAdd;
  std::function<void()> callback = nullptr;
};

QQt(std::vector<hlong> ids, const timingConfig& config, oogs_mode gsMode, const MPI_Comm& comm);
QQt(oogs_t *h);

void startFinish(const std::string& op_,
                 occa::memory& o_v,
                 const dlong stride,
                 const int k);

void startFinish(const std::string& op_,
                 occa::memory& o_v,
                 const dlong stride);

void startFinish(const std::string& op_,
                 occa::memory& o_v);

private:
  oogs_t *gsH;

};

#endif
