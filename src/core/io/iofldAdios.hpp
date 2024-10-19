#ifdef NEKRS_ENABLE_ADIOS

#if !defined(nekrs_iofldAdios_hpp_)
#define nekrs_iofldAdios_hpp_

#include "iofld.hpp"
#include "adios2.h"
#include "nekInterfaceAdapter.hpp"
#include "pointInterpolation.hpp"

class iofldAdios : public iofld
{
public:
  ~iofldAdios()
  {
    close();
  };

  void openEngine() override;
  size_t write() override;
  size_t read() override;
  void close() override;

private:
  static constexpr int VTK_HEXAHEDRON = 12;
  static constexpr const char *configFile = "adios.yaml";

  std::unique_ptr<pointInterpolation_t> interp;

  std::string streamName;

  adios2::ADIOS *adios;
  adios2::IO adiosIO;
  adios2::Engine adiosEngine;

  template <typename OutputType> size_t write_();

  template <typename InputType, typename OutputType>
  void putVariableConvert(const std::vector<occa::memory> &o_fld, occa::memory &fldConv);

  void generateConnectivity(occa::memory etov);
  std::string vtkSchema();

  void validateUserSingleValues(const std::string &name) override;
  void validateUserFields(const std::string &name) override;

  template <typename T>
  void putVariable(const std::string &name,
                   const occa::memory &data,
                   adios2::Dims dim,
                   adios2::Mode putMode = adios2::Mode::Deferred)
  {
    auto var = (adiosIO.InquireVariable<T>(name))
                   ? adiosIO.InquireVariable<T>(name)
                   : adiosIO.DefineVariable<T>(name, {}, {}, dim, adios2::ConstantDims);

    adiosEngine.Put(var, static_cast<const T *>(data.ptr()), putMode);
  };

  template <typename T>
  void putVariable(const std::string &name, const T &data, adios2::Mode putMode = adios2::Mode::Deferred)
  {
    auto var = (adiosIO.InquireVariable<T>(name)) ? adiosIO.InquireVariable<T>(name)
                                                  : adiosIO.DefineVariable<T>(name);

    adiosEngine.Put(var, data, putMode);
  };

  void
  putVariable(std::string name, const variantType &variant, adios2::Mode putMode = adios2::Mode::Deferred)
  {
    std::visit(
        [this, &name](const auto &variant) {
          if constexpr (std::is_same_v<std::decay_t<decltype(variant)>, std::reference_wrapper<int>>) {
            putVariable<int>(name, variant.get());
          } else if constexpr (std::is_same_v<std::decay_t<decltype(variant)>,
                                              std::reference_wrapper<long long int>>) {
            putVariable<long long int>(name, variant.get());
          } else if constexpr (std::is_same_v<std::decay_t<decltype(variant)>,
                                              std::reference_wrapper<float>>) {
            putVariable<float>(name, variant.get());
          } else if constexpr (std::is_same_v<std::decay_t<decltype(variant)>,
                                              std::reference_wrapper<double>>) {
            putVariable<double>(name, variant.get());
          }
        },
        variant);
  };

  typedef struct {
    std::string type;
    int dim;
    size_t step;
    std::vector<std::pair<size_t, size_t>> blocks;
    occa::memory data;
  } variable;

  std::map<std::string, variable> variables;

  template <class T> int getVariable(bool registerOnly, const std::string &name, size_t step);

  template <typename Tin, typename Tout> std::vector<occa::memory> getDataConvert(const std::string &name);

  template <typename T> std::vector<occa::memory> redistributeField(const std::vector<occa::memory> &o_in);

  template <typename T> void getData(const std::string &name, std::vector<occa::memory> &o_out);

  template <typename T> void getData(const std::string &name, variantType &variant);
};

#endif

#endif
