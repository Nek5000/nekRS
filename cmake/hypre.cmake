set(HYPRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/hypre)

set(HYPRE_FLAGS_EXTRA "-fPIC")
if("${CMAKE_CXX_COMPILE_OPTIONS_VISIBILITY}" STREQUAL "")
  if(USING_NVHPC)
    # CMake doesn't populate this flag correctly - at least for the moment
    set(CMAKE_CXX_COMPILE_OPTIONS_VISIBILITY "-fvisibility=")
  else()
    message(WARNING "Cannot identify visibility compiler flag!")
  endif()
endif()
string(APPEND HYPRE_FLAGS_EXTRA " ${CMAKE_CXX_COMPILE_OPTIONS_VISIBILITY}hidden")


set(HYPRE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/HYPRE_BUILD-prefix)
set(HYPRE_BUILD_DIR ${HYPRE_INSTALL_DIR}/src/HYPRE_BUILD)
ExternalProject_Add(
   HYPRE_BUILD
   URL "${HYPRE_SOURCE_DIR}" 
   CONFIGURE_COMMAND cd ${HYPRE_BUILD_DIR}/src && ./configure 
     --prefix=${HYPRE_INSTALL_DIR}
     --with-extra-CFLAGS=${HYPRE_FLAGS_EXTRA}
     --with-extra-CXXFLAGS=${HYPRE_FLAGS_EXTRA}
     --disable-shared --enable-single --enable-mixedint --disable-fortran
     ${HYPRE_CONFIGURE_FLAGS}
   BUILD_COMMAND "" 
   INSTALL_COMMAND cd ${HYPRE_BUILD_DIR}/src && $(MAKE) install
)

add_library(nekrs-hypre SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/hypre/hypreWrapper.cpp)
add_dependencies(nekrs-hypre HYPRE_BUILD)
target_include_directories(nekrs-hypre PRIVATE ${HYPRE_INSTALL_DIR}/include)
target_link_libraries(nekrs-hypre PUBLIC MPI::MPI_C 
                                  PRIVATE ${HYPRE_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}HYPRE.a)
set_target_properties(nekrs-hypre PROPERTIES CXX_VISIBILITY_PRESET hidden)


if(ENABLE_HYPRE_GPU AND (OCCA_CUDA_ENABLED OR OCCA_HIP_ENABLED))

if(OCCA_CUDA_ENABLED)
  enable_language(CUDA)
  find_package(CUDAToolkit 11.0 REQUIRED)

  set(HYPRE_DEVICE_COMPILER "${CMAKE_CUDA_COMPILER} -ccbin=${CMAKE_CXX_COMPILER}")
  set(HYPRE_COMPILER_C_FLAGS ${HYPRE_FLAGS_EXTRA})
  set(HYPRE_COMPILER_CXX_FLAGS ${HYPRE_FLAGS_EXTRA})
  set(HYPRE_DEVICE_COMPILER_FLAGS "")

  set(HYPRE_DEP "CUDA::cudart")
  list(APPEND HYPRE_DEP "CUDA::curand")
  list(APPEND HYPRE_DEP "CUDA::cublas")
  list(APPEND HYPRE_DEP "CUDA::cusparse")
  list(APPEND HYPRE_DEP "CUDA::cusolver")

  set(HYPRE_BACKEND "--with-cuda" "--with-cuda-home=${CUDAToolkit_LIBRARY_ROOT}")

  if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL "12.0.0")
    set(HYPRE_DEVICE_ARCH "HYPRE_CUDA_SM=80 90")
    #disable for now as it might not play well with all MPI implementations
    #set(HYPRE_CONFIGURE_FLAGS "--enable-device-malloc-async")
  else()
    set(HYPRE_DEVICE_ARCH "HYPRE_CUDA_SM=70 80")
  endif()

elseif(OCCA_HIP_ENABLED)
  enable_language(HIP)

  set(HYPRE_DEVICE_COMPILER ${CMAKE_CXX_COMPILER}) # used for device + host code
  set(HYPRE_COMPILER_C_FLAGS ${HYPRE_FLAGS_EXTRA})
  set(HYPRE_COMPILER_CXX_FLAGS ${HYPRE_FLAGS_EXTRA})
  set(HYPRE_DEVICE_COMPILER_FLAGS ${HYPRE_FLAGS_EXTRA})

  find_package(rocrand REQUIRED)
  find_package(rocblas REQUIRED)
  find_package(rocsparse REQUIRED)
  find_package(rocsolver REQUIRED)

  set(HYPRE_DEP "roc::rocrand")
  list(APPEND HYPRE_DEP "roc::rocblas")
  list(APPEND HYPRE_DEP "roc::rocsparse")
  list(APPEND HYPRE_DEP "roc::rocsolver")

  set(HYPRE_BACKEND "--with-hip")
endif()

if(NEKRS_GPU_MPI)
  list(APPEND HYPRE_CONFIGURE_FLAGS --enable-gpu-aware-mpi)
endif()

  set(HYPRE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/HYPRE_BUILD_DEVICE-prefix)
  set(HYPRE_BUILD_DIR ${HYPRE_INSTALL_DIR}/src/HYPRE_BUILD_DEVICE)
  ExternalProject_Add(
   HYPRE_BUILD_DEVICE
   URL "${HYPRE_SOURCE_DIR}" 
   CONFIGURE_COMMAND cd ${HYPRE_BUILD_DIR}/src && ./configure
     CUCC=${HYPRE_DEVICE_COMPILER}
     --prefix=${HYPRE_INSTALL_DIR}
     --with-extra-CFLAGS=${HYPRE_COMPILER_C_FLAGS}
     --with-extra-CXXFLAGS=${HYPRE_COMPILER_CXX_FLAGS}
     --with-extra-CUFLAGS=${HYPRE_DEVICE_COMPILER_FLAGS}
     --disable-shared --enable-single --enable-mixedint --disable-fortran
     ${HYPRE_BACKEND} ${HYPRE_DEVICE_ARCH}
     ${HYPRE_CONFIGURE_FLAGS}
   BUILD_COMMAND "" 
   INSTALL_COMMAND cd ${HYPRE_BUILD_DIR}/src && $(MAKE) install
  )

  add_library(nekrs-hypre-device SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/hypre/hypreWrapperDevice.cpp)
  add_dependencies(nekrs-hypre-device HYPRE_BUILD_DEVICE)
  target_compile_definitions(nekrs-hypre-device PRIVATE -DENABLE_HYPRE_GPU)
  target_include_directories(nekrs-hypre-device PRIVATE ${HYPRE_INSTALL_DIR}/include)
  target_link_libraries(nekrs-hypre-device 
                        PUBLIC libocca MPI::MPI_C 
                        PRIVATE ${HYPRE_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}HYPRE.a 
			${HYPRE_DEP}) 
  set_target_properties(nekrs-hypre-device PROPERTIES CXX_VISIBILITY_PRESET hidden)

else()
  #dummy
  message(WARNING "No supported HYPRE backend found - disable device support!")
  add_library(nekrs-hypre-device SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/hypre/hypreWrapperDevice.cpp)
  target_link_libraries(nekrs-hypre-device PUBLIC libocca MPI::MPI_C) 
endif()

unset(HYPRE_DEVICE_COMPILER)
unset(HYPRE_BUILD_DIR)
unset(HYPRE_SOURCE_DIR)
unset(HYPRE_CONFIGURE_FLAGS) 
unset(HYPRE_INSTALL_DIR)
unset(HYPRE_FLAGS_EXTRA)
unset(HYPRE_COMPILER_C_FLAGS)
unset(HYPRE_COMPILER_CXX_FLAGS)
unset(HYPRE_DEVICE_COMPILER_FLAGS)
unset(HYPRE_DEP)
