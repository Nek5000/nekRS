set(HYPRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/hypre)

set(HYPRE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/HYPRE_BUILD-prefix)
ExternalProject_Add(
    HYPRE_BUILD
    SOURCE_DIR ${HYPRE_SOURCE_DIR}
    SOURCE_SUBDIR "src"
#    BUILD_ALWAYS ON
    CMAKE_ARGS  -DHYPRE_INSTALL_PREFIX=${HYPRE_INSTALL_DIR}
                -DHYPRE_BUILD_TYPE=RelWithDebInfo
                -DCMAKE_INSTALL_LIBDIR=${HYPRE_INSTALL_DIR}/lib
                -DCMAKE_C_FLAGS_RELWITHDEBINFO=${CMAKE_C_FLAGS_RELWITHDEBINFO}
                -DHYPRE_ENABLE_SHARED=OFF
                -DHYPRE_ENABLE_MIXEDINT=ON
                -DHYPRE_ENABLE_SINGLE=ON
                -DHYPRE_WITH_OPENMP=OFF
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                -DCMAKE_C_VISIBILITY_PRESET=hidden
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
)

add_library(nekrs-hypre SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/hypre/hypreWrapper.cpp)
add_dependencies(nekrs-hypre HYPRE_BUILD)
target_include_directories(nekrs-hypre PRIVATE ${HYPRE_INSTALL_DIR}/include)
# lacking of a better alternative adding dependencies manually 
target_link_libraries(nekrs-hypre PUBLIC MPI::MPI_C 
                                  PRIVATE ${HYPRE_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}HYPRE.a)
set_target_properties(nekrs-hypre PROPERTIES CXX_VISIBILITY_PRESET hidden)

if(ENABLE_HYPRE_GPU)

if(ENABLE_CUDA)
  find_package(CUDAToolkit 10.0 REQUIRED)
  set(HYPRE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/HYPRE_BUILD_DEVICE-prefix)
  set(HYPRE_CUDA_SM 70)

  if(CUDA_VERSION VERSION_GREATER_EQUAL 11.1.0)
    set(HYPRE_CUDA_SM 70 80)
  endif()

#  if(CUDA_VERSION VERSION_GREATER_EQUAL 11.2.0)
#    set(HYPRE_ENABLE_DEVICE_MALLOC_ASYNC ON)
#  endif()

  ExternalProject_Add(
      HYPRE_BUILD_DEVICE
      SOURCE_DIR ${HYPRE_SOURCE_DIR}
      SOURCE_SUBDIR "src"
#      BUILD_ALWAYS ON
      CMAKE_CACHE_ARGS
                  -DHYPRE_CUDA_SM:STRING=${HYPRE_CUDA_SM}
      CMAKE_ARGS  
                  -DHYPRE_ENABLE_SHARED=OFF
                  -DHYPRE_ENABLE_MIXEDINT=ON
                  -DHYPRE_ENABLE_SINGLE=ON
                  -DHYPRE_WITH_CUDA=ON
                  -DHYPRE_WITH_GPU_AWARE_MPI=${NEKRS_GPU_MPI}
                  -DHYPRE_ENABLE_CUSPARSE=ON
                  -DHYPRE_ENABLE_DEVICE_MALLOC_ASYNC=${HYPRE_ENABLE_DEVICE_MALLOC_ASYNC}
                  -DHYPRE_BUILD_TYPE=RelWithDebInfo
                  -DHYPRE_INSTALL_PREFIX=${HYPRE_INSTALL_DIR}
                  -DCMAKE_INSTALL_LIBDIR=${HYPRE_INSTALL_DIR}/lib
                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  -DCMAKE_C_FLAGS_RELWITHDEBINFO=${CMAKE_C_FLAGS_RELWITHDEBINFO}
                  -DCMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
                  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                  -DCMAKE_C_VISIBILITY_PRESET=hidden
                  -DCMAKE_CXX_VISIBILITY_PRESET=hidden
                  -DCMAKE_CUDA_VISIBILITY_PRESET=hidden
                  -DCMAKE_CUDA_HOST_COMPILER=${CMAKE_CXX_COMPILER}

  )

  add_library(nekrs-hypre-device SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/hypre/hypreWrapperDevice.cpp)
  add_dependencies(nekrs-hypre-device HYPRE_BUILD_DEVICE)
  target_compile_definitions(nekrs-hypre-device PRIVATE -DENABLE_HYPRE_GPU)
  target_include_directories(nekrs-hypre-device PRIVATE ${HYPRE_INSTALL_DIR}/include)
  # lacking of a better alternative adding dependencies manually 
  target_link_libraries(nekrs-hypre-device 
                        PUBLIC libocca MPI::MPI_C 
                        PRIVATE ${HYPRE_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}HYPRE.a 
                        CUDA::curand CUDA::cublas CUDA::cusparse CUDA::cusolver) 
  set_target_properties(nekrs-hypre-device PROPERTIES CXX_VISIBILITY_PRESET hidden)
elseif(ENABLE_HIP)
  message(FATAL_ERROR "HYPRE wrapper build does not support HIP!")
endif()

else()
  #dummy
  add_library(nekrs-hypre-device SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/elliptic/amgSolver/hypre/hypreWrapperDevice.cpp)
  target_link_libraries(nekrs-hypre-device PUBLIC libocca MPI::MPI_C) 
endif()
