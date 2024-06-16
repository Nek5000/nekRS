#include <stdio.h>

#include <occa/internal/modes/dpcpp/utils.hpp>
#include <occa/internal/modes/dpcpp/device.hpp>
#include <occa/internal/modes/dpcpp/memory.hpp>
#include <occa/internal/modes/dpcpp/kernel.hpp>
#include <occa/internal/modes/dpcpp/stream.hpp>
#include <occa/internal/modes/dpcpp/streamTag.hpp>
#include <occa/internal/io.hpp>
#include <occa/internal/utils/sys.hpp>
#include <occa/internal/utils/env.hpp>
#include <occa/core/base.hpp>

namespace occa
{
  namespace dpcpp
  {
    /* Returns true if any DPC++ device is enabled on the machine */
    bool isEnabled()
    {
      auto device_list = ::sycl::device::get_devices();
      return (device_list.size() > 0);
    }

    void setCompiler(occa::json &dpcpp_properties) noexcept
    {
      std::string compiler;
      if (env::var("OCCA_DPCPP_COMPILER").size()) {
        compiler = env::var("OCCA_DPCPP_COMPILER");
      } else if (dpcpp_properties.has("compiler")){
        compiler = dpcpp_properties["compiler"].toString();
      } else if (env::var("OCCA_CXX").size()) {
        compiler = env::var("OCCA_CXX");
      } else if (env::var("CXX").size()) {
        compiler = env::var("CXX");
      } else {
        OCCA_FORCE_WARNING("OCCA_DPCPP_COMPILER is defaulting to clang++");
        compiler = "clang++";
      }
      dpcpp_properties["compiler"] = compiler;
    }

    void setCompilerFlags(occa::json &dpcpp_properties) noexcept
    {
      std::string compiler_flags;
      if (dpcpp_properties.has("compiler_flags"))
      {
        compiler_flags = dpcpp_properties["compiler_flags"].toString();
      }
      else if (env::var("OCCA_DPCPP_COMPILER_FLAGS").size())
      {
        compiler_flags = env::var("OCCA_DPCPP_COMPILER_FLAGS");
      }
      else
      {
        compiler_flags = "-O3 -fsycl";
      }
      dpcpp_properties["compiler_flags"] = compiler_flags;
    }

    void setSharedFlags(occa::json &dpcpp_properties) noexcept
    {
      std::string shared_flags;
      if (env::var("OCCA_COMPILER_SHARED_FLAGS").size())
      {
        shared_flags = env::var("OCCA_COMPILER_SHARED_FLAGS");
      }
      else if (dpcpp_properties.has("compiler_shared_flags"))
      {
        shared_flags = (std::string) dpcpp_properties["compiler_shared_flags"];
      }
      else
      {
        shared_flags = "-shared -fPIC";
      }
      dpcpp_properties["compiler_shared_flags"] = shared_flags;
    }

    void setLinkerFlags(occa::json &dpcpp_properties) noexcept
    {
      std::string linker_flags;
      if (env::var("OCCA_DPCPP_LINKER_FLAGS").size())
      {
        linker_flags = env::var("OCCA_DPCPP_LINKER_FLAGS");
      }
      else if (dpcpp_properties.has("linker_flags"))
      {
        linker_flags = dpcpp_properties["linker_flags"].toString();
      }
      dpcpp_properties["linker_flags"] = linker_flags;
    }

    std::string getDeviceArchString(const ::sycl::device& sycl_device) {
    #if SYCL_EXT_ONEAPI_DEVICE_ARCHITECTURE
      namespace syclex = sycl::ext::oneapi::experimental;
      using arch_type = syclex::architecture;
      arch_type arch = sycl_device.get_info<syclex::info::device::architecture>();
      switch (arch) {
        case arch_type::intel_gpu_pvc     : return "intel_gpu_pvc";
        case arch_type::intel_gpu_acm_g12 : return "intel_gpu_acm_g12";
        case arch_type::intel_gpu_acm_g11 : return "intel_gpu_acm_g11";
        case arch_type::intel_gpu_acm_g10 : return "intel_gpu_acm_g10";
        case arch_type::intel_gpu_adl_n   : return "intel_gpu_adl_n";
        case arch_type::intel_gpu_adl_p   : return "intel_gpu_adl_p";
        case arch_type::intel_gpu_adl_s   : return "intel_gpu_adl_s";
        case arch_type::intel_gpu_rkl     : return "intel_gpu_rkl";
        case arch_type::intel_gpu_ehl     : return "intel_gpu_ehl";
        case arch_type::intel_gpu_12_10_0 : return "intel_gpu_12_10";
        case arch_type::intel_gpu_12_0_0  : return "intel_gpu_12_0_";
        case arch_type::intel_gpu_11_0_0  : return "intel_gpu_11_0_0";
        case arch_type::intel_gpu_9_7_0   : return "intel_gpu_9_7_0";
        case arch_type::intel_gpu_9_6_0   : return "intel_gpu_9_6_0";
        case arch_type::intel_gpu_9_5_0   : return "intel_gpu_9_5_0";
        case arch_type::intel_gpu_9_4_0   : return "intel_gpu_9_4_0";
        case arch_type::intel_gpu_9_3_0   : return "intel_gpu_9_3_0";
        case arch_type::intel_gpu_9_2_9   : return "intel_gpu_9_2_9";
        case arch_type::intel_gpu_9_1_9   : return "intel_gpu_9_1_9";
        case arch_type::intel_gpu_9_0_9   : return "intel_gpu_9_0_9";
        case arch_type::intel_gpu_8_0_0   : return "intel_gpu_8_0_0";
        case arch_type::nvidia_gpu_sm_50  : return "nvidia_gpu_sm_50";
        case arch_type::nvidia_gpu_sm_52  : return "nvidia_gpu_sm_52";
        case arch_type::nvidia_gpu_sm_53  : return "nvidia_gpu_sm_53";
        case arch_type::nvidia_gpu_sm_60  : return "nvidia_gpu_sm_60";
        case arch_type::nvidia_gpu_sm_61  : return "nvidia_gpu_sm_61";
        case arch_type::nvidia_gpu_sm_62  : return "nvidia_gpu_sm_62";
        case arch_type::nvidia_gpu_sm_70  : return "nvidia_gpu_sm_70";
        case arch_type::nvidia_gpu_sm_72  : return "nvidia_gpu_sm_72";
        case arch_type::nvidia_gpu_sm_75  : return "nvidia_gpu_sm_75";
        case arch_type::nvidia_gpu_sm_80  : return "nvidia_gpu_sm_80";
        case arch_type::nvidia_gpu_sm_86  : return "nvidia_gpu_sm_86";
        case arch_type::nvidia_gpu_sm_87  : return "nvidia_gpu_sm_87";
        case arch_type::nvidia_gpu_sm_89  : return "nvidia_gpu_sm_89";
        case arch_type::nvidia_gpu_sm_90  : return "nvidia_gpu_sm_90";   
        case arch_type::amd_gpu_gfx700    : return "amd_gpu_gfx700";
        case arch_type::amd_gpu_gfx701    : return "amd_gpu_gfx701";
        case arch_type::amd_gpu_gfx702    : return "amd_gpu_gfx702";
        case arch_type::amd_gpu_gfx801    : return "amd_gpu_gfx801";
        case arch_type::amd_gpu_gfx802    : return "amd_gpu_gfx802";
        case arch_type::amd_gpu_gfx803    : return "amd_gpu_gfx803";
        case arch_type::amd_gpu_gfx805    : return "amd_gpu_gfx805";
        case arch_type::amd_gpu_gfx810    : return "amd_gpu_gfx810";
        case arch_type::amd_gpu_gfx900    : return "amd_gpu_gfx900"; 
        case arch_type::amd_gpu_gfx902    : return "amd_gpu_gfx902"; 
        case arch_type::amd_gpu_gfx904    : return "amd_gpu_gfx904"; 
        case arch_type::amd_gpu_gfx906    : return "amd_gpu_gfx906"; 
        case arch_type::amd_gpu_gfx908    : return "amd_gpu_gfx908"; 
        // case arch_type::amd_gpu_gfx909    : return "amd_gpu_gfx909"; 
        case arch_type::amd_gpu_gfx90a    : return "amd_gpu_gfx90a"; 
        // case arch_type::amd_gpu_gfx90c    : return "amd_gpu_gfx90c"; 
        // case arch_type::amd_gpu_gfx940    : return "amd_gpu_gfx940"; 
        // case arch_type::amd_gpu_gfx941    : return "amd_gpu_gfx941"; 
        // case arch_type::amd_gpu_gfx942    : return "amd_gpu_gfx942"; 
        case arch_type::amd_gpu_gfx1010   : return "amd_gpu_gfx1010";
        case arch_type::amd_gpu_gfx1011   : return "amd_gpu_gfx1011";
        case arch_type::amd_gpu_gfx1012   : return "amd_gpu_gfx1012";
        case arch_type::amd_gpu_gfx1013   : return "amd_gpu_gfx1013";
        case arch_type::amd_gpu_gfx1030   : return "amd_gpu_gfx1030";
        case arch_type::amd_gpu_gfx1031   : return "amd_gpu_gfx1031";
        case arch_type::amd_gpu_gfx1032   : return "amd_gpu_gfx1032";
        // case arch_type::amd_gpu_gfx1033   : return "amd_gpu_gfx1033";
        case arch_type::amd_gpu_gfx1034   : return "amd_gpu_gfx1034";
        // case arch_type::amd_gpu_gfx1035   : return "amd_gpu_gfx1035";
        // case arch_type::amd_gpu_gfx1036   : return "amd_gpu_gfx1036";
        // case arch_type::amd_gpu_gfx1100   : return "amd_gpu_gfx1100"; 
        // case arch_type::amd_gpu_gfx1101   : return "amd_gpu_gfx1101"; 
        // case arch_type::amd_gpu_gfx1102   : return "amd_gpu_gfx1102"; 
        // case arch_type::amd_gpu_gfx1103   : return "amd_gpu_gfx1103"; 
        // case arch_type::amd_gpu_gfx1150   : return "amd_gpu_gfx1150"; 
        // case arch_type::amd_gpu_gfx1151   : return "amd_gpu_gfx1151"; 
        // case arch_type::amd_gpu_gfx1200   : return "amd_gpu_gfx1200"; 
        // case arch_type::amd_gpu_gfx1201   : return "amd_gpu_gfx1201"; 
        default: return "";
      }
    #else
      return "";
    #endif
    }

    void setDeviceArchCompilerFlags(const ::sycl::device& sycl_device, 
                                    std::string& compiler_flags) {

      const std::string arch_string = getDeviceArchString(sycl_device);
      if (std::string::npos == compiler_flags.find("-fsycl-targets")) {
        if (!arch_string.empty()) {
          compiler_flags += " -fsycl-targets=" + arch_string;
        } 
      }
      if (std::string::npos == compiler_flags.find("-ftarget-register-alloc-mode")) {
        if (std::string::npos != arch_string.find("pvc")) {
          compiler_flags += " -Xsycl-target-backend '-options -ze-intel-enable-auto-large-GRF-mode'";
        }
      }
    }

    occa::dpcpp::device& getDpcppDevice(modeDevice_t* device_)
    {
      occa::dpcpp::device* dpcppDevice = dynamic_cast<occa::dpcpp::device*>(device_);
      OCCA_ERROR("[dpcpp::getDpcppDevice] Dynamic cast failed!",nullptr != dpcppDevice);
      return *dpcppDevice;
    }  

    occa::dpcpp::stream& getDpcppStream(const occa::stream& stream_)
    {
      auto* dpcpp_stream{dynamic_cast<occa::dpcpp::stream*>(stream_.getModeStream())};
      OCCA_ERROR("[dpcpp::getDpcppStream]: Dynamic cast failed!", nullptr != dpcpp_stream);
      return *dpcpp_stream;
    }

    occa::dpcpp::streamTag& getDpcppStreamTag(const occa::streamTag& tag_)
    {
      auto* dpcppTag{dynamic_cast<occa::dpcpp::streamTag*>(tag_.getModeStreamTag())};
      OCCA_ERROR("[dpcpp::getDpcppStreamTag]: Dynamic cast failed!", nullptr != dpcppTag);
      return *dpcppTag;
    }

    occa::device wrapDevice(::sycl::device sycl_device,
                            const occa::properties &props)
    {
      occa::properties allProps;
      allProps["mode"] = "dpcpp";
      allProps["wrapped"] = true;
      allProps += props;

      auto* wrapper{new dpcpp::device(allProps, sycl_device)};
      wrapper->dontUseRefs();

      wrapper->currentStream = wrapper->createStream(allProps["stream"]);
      return occa::device(wrapper);
    }

    void warn(const ::sycl::exception &e,
              const std::string &filename,
              const std::string &function,
              const int line,
              const std::string &message)
    {
      std::stringstream ss;
      ss << message << "\n"
         << "DPCPP Error:"
         << e.what();
      occa::warn(filename, function, line, ss.str());
    }

    void error(const ::sycl::exception &e,
               const std::string &filename,
               const std::string &function,
               const int line,
               const std::string &message)
    {
      std::stringstream ss;
      ss << message << "\n"
         << "DPCPP Error:"
         << e.what();
      occa::error(filename, function, line, ss.str());
    }

  } // namespace dpcpp
} // namespace occa
