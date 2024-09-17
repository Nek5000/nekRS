#if !defined(nekrs_iofldNek_hpp_)
#define nekrs_iofldNek_hpp_

#include "iofld.hpp"
#include "nekInterfaceAdapter.hpp"

class iofldNek : public iofld
{
public:
  ~iofldNek()
  {
    close();
  };

  void validateUserFields(const std::string &name) override;
  void validateUserSingleValues(const std::string &name) override;
  void openEngine() override;
  size_t write() override;
  size_t read() override;
  void close() override;

private:
  nek::fldData fldData;
  std::string fileSuffix();
};

#endif
