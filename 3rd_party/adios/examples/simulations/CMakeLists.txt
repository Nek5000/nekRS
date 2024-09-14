if(ADIOS2_HAVE_MPI)
  add_subdirectory(gray-scott)
  add_subdirectory(gray-scott-struct)

  if(ADIOS2_HAVE_Kokkos)
    add_subdirectory(gray-scott-kokkos)
  endif()

  add_subdirectory(heatTransfer)
endif()

add_subdirectory(korteweg-de-vries)
add_subdirectory(lorenz_ode)
