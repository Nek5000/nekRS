cmake_minimum_required(VERSION 3.22)
project(gemv VERSION 0.0.1 DESCRIPTION "Platform independent gemv on GPUs"
  LANGUAGES C)

# Build Options:
option(ENABLE_CUDA "Build CUDA backend" ON)
option(ENABLE_HIP "Build HIP backend" ON)

# Build defaults (which can be overridden at runtime.
set(GEMV_DEFAULT_VERBOSE 1)
set(GEMV_DEFAULT_DEVICE 0)
set(GEMV_DEFAULT_PRECISION 0)
set(GEMV_DEFAULT_BACKEND "")
set(GEMV_MAX_BACKEND_LENGTH 16)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(GEMV_SOURCES src/gemv.c src/gemv-impl.c src/gemv-log.c
  backends/gemv-backend-dummy.c)

if (ENABLE_HIP)
  find_package(HIP)
  if (TARGET gemv::HIP)
    set(GEMV_SOURCES ${GEMV_SOURCES} backends/gemv-backend-hip.c)
  endif()
  if (TARGET roc::hipblas)
    set(GEMV_SOURCES ${GEMV_SOURCES} backends/gemv-backend-hipblas.c)
  endif()
endif()

if (ENABLE_CUDA)
  find_package(CUDA)
  if (TARGET gemv::CUDA)
    set(GEMV_SOURCES ${GEMV_SOURCES} backends/gemv-backend-cuda.c)
  endif()
  if (TARGET cuda::cublas)
    set(GEMV_SOURCES ${GEMV_SOURCES} backends/gemv-backend-cublas.c)
  endif()
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/gemv-defs.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/include/gemv-defs.h @ONLY)

add_library(gemv SHARED ${GEMV_SOURCES})
set_target_properties(gemv PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION 1
  PUBLIC_HEADER include/gemv.h)
target_include_directories(gemv PRIVATE include)
target_compile_options(gemv PRIVATE $<$<C_COMPILER_ID:MSVC>:/W4 /WX>
  $<$<NOT:$<C_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic -Wno-unused-function
  -Wno-unused-parameter -Wno-stringop-overflow>)

if (ENABLE_HIP AND TARGET gemv::HIP)
  target_link_libraries(gemv PRIVATE gemv::HIP)
endif()

if (ENABLE_CUDA AND TARGET gemv::CUDA)
  target_link_libraries(gemv PRIVATE gemv::CUDA)
endif()

macro(gemv_add_binary file_name)
  get_filename_component(target ${file_name} NAME_WE)
  set(target ${target}-bin)
  add_executable(${target} ${file_name})
  set_target_properties(${target} PROPERTIES VERSION ${PROJECT_VERSION})
  target_compile_options(${target} PRIVATE $<$<C_COMPILER_ID:MSVC>:/W4 /MX>
    $<$<NOT:$<C_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic -Wno-unused-function
    -Wno-unused-parameter -Wno-stringop-overflow>)
  target_include_directories(${target} PRIVATE include)
  target_link_libraries(${target} PRIVATE gemv)
  install(TARGETS ${target} DESTINATION bin)
endmacro()

gemv_add_binary(bin/gemv.c)

install(TARGETS gemv LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include)

# Add clang-format as a custom target if available.
find_program(CLANG_FORMAT NAMES clang-format)
if (CLANG_FORMAT)
  set(FILES_TO_FORMAT ${CMAKE_SOURCE_DIR}/src/*.c ${CMAKE_SOURCE_DIR}/backends/*.c ${CMAKE_SOURCE_DIR}/bin/*.c)
  set(FILES_TO_FORMAT ${FILES_TO_FORMAT} ${CMAKE_SOURCE_DIR}/include/*.h ${CMAKE_SOURCE_DIR}/include/backends/*.h)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT} -i ${FILES_TO_FORMAT}
    COMMENT "Running clang-format ...")
  add_custom_target(format-check
    COMMAND ${CLANG_FORMAT} --dry-run -Werror -i ${FILES_TO_FORMAT}
    COMMENT "Running clang-format check ...")
endif()

# Configure .clang-tidy.in to generate .clang-tidy.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy.in
  ${CMAKE_CURRENT_SOURCE_DIR}/.clang-tidy @ONLY)
# Add clang-tidy as a custom target if available.
find_program(CLANG_TIDY NAMES clang-tidy)
if (CLANG_TIDY)
  add_custom_target(tidy
    COMMAND ${CLANG_TIDY} -p ${CMAKE_CURRENT_BINARY_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/*/*.[ch]
    COMMENT "Running clang-tidy ...")
endif()
