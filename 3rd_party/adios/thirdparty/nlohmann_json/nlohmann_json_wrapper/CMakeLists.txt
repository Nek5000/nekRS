add_library(nlohmann_json_wrapper INTERFACE)
target_sources(nlohmann_json_wrapper INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nlohmann_json.hpp>
)
target_include_directories(nlohmann_json_wrapper INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_compile_options(nlohmann_json_wrapper INTERFACE
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Intel>>:-Wno-all -Wno-error>
)

if(ADIOS2_USE_EXTERNAL_NLOHMANN_JSON)
  find_package(nlohmann_json REQUIRED)
  target_link_libraries(nlohmann_json_wrapper INTERFACE
    $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>
  )
else()
  target_sources(nlohmann_json_wrapper INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/single_include/nlohmann/json.hpp>
  )
  target_include_directories(nlohmann_json_wrapper INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/single_include>
  )
endif()

install(TARGETS nlohmann_json_wrapper EXPORT adios2Exports)
