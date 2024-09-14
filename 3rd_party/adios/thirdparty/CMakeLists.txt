# Disable warnings in thirdparty code
if(CMAKE_C_COMPILER_ID MATCHES
    "^(GNU|Clang|AppleClang|XLClang|XL|VisualAge|SunPro|HP|Intel|IntelLLVM|PGI|NVIDIA|NVHPC)$")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES
    "^(GNU|Clang|AppleClang|XLClang|XL|VisualAge|SunPro|HP|Intel|IntelLLVM|PGI|NVIDIA|NVHPC)$")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()


function(message_start_thirdparty)
  get_filename_component(tp_base "${CMAKE_CURRENT_LIST_DIR}" NAME_WE)
  message("")
  message(STATUS "ADIOS2 ThirdParty: Configuring ${tp_base}")
endfunction()

function(message_end_thirdparty)
endfunction()

include(CMakeDependentOption)

option(ADIOS2_USE_EXTERNAL_DEPENDENCIES
  "Use externally supplied dependencies" OFF)
option(ADIOS2_USE_EXTERNAL_GTEST
  "Use an externally supplied GTest library" "${ADIOS2_USE_EXTERNAL_DEPENDENCIES}")
option(ADIOS2_USE_EXTERNAL_PUGIXML
  "Use an externally supplied pugixml library" "${ADIOS2_USE_EXTERNAL_DEPENDENCIES}")
option(ADIOS2_USE_EXTERNAL_PYBIND11
  "Use an externally supplied pybind11 library" "${ADIOS2_USE_EXTERNAL_DEPENDENCIES}")
option(ADIOS2_USE_EXTERNAL_NLOHMANN_JSON
  "Use an externally supplied nlohmann_json library" "${ADIOS2_USE_EXTERNAL_DEPENDENCIES}")
option(ADIOS2_USE_EXTERNAL_YAMLCPP
  "Use an externally supplied yaml-cpp library" "${ADIOS2_USE_EXTERNAL_DEPENDENCIES}")

add_subdirectory(KWSys)

add_subdirectory(perfstubs)
adios2_add_thirdparty_target(perfstubs-interface perfstubs-interface)

if(BUILD_TESTING)
  if(ADIOS2_USE_EXTERNAL_GTEST)
    find_package(GTest REQUIRED)
    adios2_add_thirdparty_target(gtest GTest::GTest)
  else()
    add_subdirectory(GTest)
    adios2_add_thirdparty_target(gtest gtest)
  endif()
endif()

if(ADIOS2_HAVE_Python)
  if(ADIOS2_USE_EXTERNAL_PYBIND11)
    find_package(pybind11 REQUIRED)
    adios2_add_thirdparty_target(pybind11 pybind11::pybind11)
  else()
    add_subdirectory(pybind11)
    adios2_add_thirdparty_target(pybind11 pybind11)
  endif()
endif()

if(ADIOS2_USE_EXTERNAL_PUGIXML)
  find_package(pugixml 1.10 REQUIRED)
else()
  add_subdirectory(pugixml)
endif()
adios2_add_thirdparty_target(pugixml pugixml)

# ADIOS uses a wrapper header for nlohmann_json to suppress noisy compiler
# warnings
add_subdirectory(nlohmann_json)
adios2_add_thirdparty_target(nlohmann_json nlohmann_json_wrapper)

if(ADIOS2_USE_EXTERNAL_YAMLCPP)
  find_package(yaml-cpp 0.7.0 REQUIRED)
else()
  add_subdirectory(yaml-cpp)
endif()

adios2_add_thirdparty_target(yaml-cpp yaml-cpp)

# YAML-cpp does not add the incdir/libdir to its target rather it exposes it
# through the cmake variables YAML_CPP_INCLUDE_DIR and YAML_CPP_LIBRARY_DIR
if(yaml-cpp_VERSION VERSION_GREATER_EQUAL 0.8)
  target_include_directories(adios2::thirdparty::yaml-cpp INTERFACE ${YAML_CPP_INCLUDE_DIR})
  target_link_directories(adios2::thirdparty::yaml-cpp INTERFACE ${YAML_CPP_LIBRARY_DIR})
endif()

if(WIN32)
  add_subdirectory(mingw-w64)
endif()

#------------------------------------------------------------------------------#
# SST / EVPath stuff
#------------------------------------------------------------------------------#
cmake_dependent_option(ADIOS2_USE_EXTERNAL_EVPATH
  "Use an externally supplied EVPath library" OFF
  "ADIOS2_HAVE_SST" OFF
)
cmake_dependent_option(ADIOS2_USE_EXTERNAL_ATL
  "Use an externally supplied ATL library" OFF
  "NOT ADIOS2_USE_EXTERNAL_EVPATH" OFF
)
cmake_dependent_option(ADIOS2_USE_EXTERNAL_FFS
  "Use an externally supplied FFS library" OFF
  "NOT ADIOS2_USE_EXTERNAL_EVPATH" OFF
)
cmake_dependent_option(ADIOS2_USE_EXTERNAL_DILL
  "Use an externally supplied DILL library" OFF
  "NOT ADIOS2_USE_EXTERNAL_FFS" OFF
)
cmake_dependent_option(ADIOS2_USE_EXTERNAL_ENET
  "Use an externally supplied ENET library" OFF
  "NOT ADIOS2_USE_EXTERNAL_EVPATH" OFF
)

if(NOT ADIOS2_USE_EXTERNAL_FFS)
  if(NOT ADIOS2_USE_EXTERNAL_ATL)
    add_subdirectory(atl)
  endif()
  find_package(atl REQUIRED)

  if(NOT ADIOS2_USE_EXTERNAL_DILL)
    add_subdirectory(dill)
  endif()
  find_package(dill REQUIRED)

  add_subdirectory(ffs)
endif()
find_package(ffs REQUIRED)
adios2_add_thirdparty_target(ffs ffs::ffs)

if(NOT ADIOS2_USE_EXTERNAL_EVPATH)
  if(NOT ADIOS2_USE_EXTERNAL_ENET)
    add_subdirectory(enet)
  endif()
  find_package(enet REQUIRED)

  add_subdirectory(EVPath)
endif()
find_package(EVPath REQUIRED)
adios2_add_thirdparty_target(EVPath EVPath::EVPath)

# Keep this at the end so we see the verbose messages at least once
set(ADIOS2_THIRDPARTY_VERBOSE_DISABLE TRUE CACHE INTERNAL "" FORCE)
message("")
