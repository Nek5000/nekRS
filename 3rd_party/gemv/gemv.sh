#!/bin/bash

function print_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --help Print this help message and exit."
  echo "  --cc <compiler> Set the compiler to use for the build."
  echo "  --cflags <compiler flags> Set the compiler flags for the build."
  echo "  --build-type <Release|Debug> Build type."
  echo "  --build-dir <build directory> Build directory."
  echo "  --install-prefix <install prefix> Install prefix."
  echo "  --install Install the project."
  echo "  --format Format the source code with clang-format."
  echo "  --format-check Check if source formatting is compliant with clang-format."
  echo "  --tidy Run clang-tidy."
}

# Set default values.
: ${GEMV_CC:=cc}
: ${GEMV_CFLAGS:="-O2 -g -mtune=native -march=native"}
: ${GEMV_BUILD_TYPE:=Release}
: ${GEMV_INSTALL_PREFIX:=`pwd`/install}
: ${GEMV_BUILD_DIR:=`pwd`/build}
: ${GEMV_INSTALL:=NO}
: ${GEMV_FORMAT:=NO}
: ${GEMV_FORMAT_CHECK:=NO}
: ${GEMV_TIDY:=NO}

# Handle command line arguments.
while [[ $# -gt 0 ]]; do
  case $1 in
    --help)
      print_help
      exit 0
      ;;
    --cc)
      GEMV_CC="$2"
      shift
      shift
      ;;
    --cflags)
      GEMV_CFLAGS="$2"
      shift
      shift
      ;;
    --build-type)
      GEMV_BUILD_TYPE="$2"
      shift
      shift
      ;;
    --build-dir)
      GEMV_BUILD_DIR="$2"
      shift
      shift
      ;;
    --install-prefix)
      GEMV_INSTALL_PREFIX="$2"
      shift
      shift
      ;;
    --install)
      GEMV_INSTALL="YES"
      shift
      ;;
    --format)
      GEMV_FORMAT="YES"
      shift
      ;;
    --format-check)
      GEMV_FORMAT_CHECK="YES"
      shift
      ;;
    --tidy)
      GEMV_TIDY="YES"
      shift
      ;;
    *)
      echo "Unknown option: $1"
      print_help
      exit 1
      ;;
  esac
done
  
mkdir -p ${GEMV_BUILD_DIR} 2> /dev/null

cmake -DCMAKE_C_COMPILER=${GEMV_CC} \
  -DCMAKE_BUILD_TYPE=${GEMV_BUILD_TYPE} \
  -DCMAKE_INSTALL_PREFIX=${GEMV_INSTALL_PREFIX} \
  -B ${GEMV_BUILD_DIR} \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -S .
  
if [[ "${GEMV_FORMAT}" == "YES" ]]; then
  cmake --build ${GEMV_BUILD_DIR} --target format -j4
fi

if [[ "${GEMV_FORMAT_CHECK}" == "YES" ]]; then
  cmake --build ${GEMV_BUILD_DIR} --target format-check -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: clang-format check failed."
    exit 1
  fi
fi

if [[ "${GEMV_TIDY}" == "YES" ]]; then
  cmake --build ${GEMV_BUILD_DIR} --target tidy -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: clang-tidy failed."
    exit 1
  fi
fi

if [[ "${GEMV_INSTALL}" == "YES" ]]; then
  cmake --build ${GEMV_BUILD_DIR} --target install -j4
  if [[ $? -ne 0 ]]; then
    echo "Error: Installing failed."
    exit 1
  fi
fi
