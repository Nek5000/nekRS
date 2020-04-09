macro(install_ src_files source install_prefix)
  list(LENGTH src_files len)
  math(EXPR counter "${len}-1")
  
  if(${len} GREATER 0)
    foreach(index RANGE ${counter})
      list(GET src_files ${index} src)
      string(REPLACE "${source}" "${install_prefix}" dest ${src})
      get_filename_component(base_dir ${dest} DIRECTORY)
      execute_process(
        COMMAND bash -c "mkdir -p ${base_dir}"
      )
      execute_process(
        COMMAND bash -c "[ ${src} -nt ${dest} ] && cp ${src} ${dest}"
      )
    endforeach()
  endif()
endmacro()

macro(install_glob_recurse_if_newer source pattern install_prefix)
  file(GLOB_RECURSE src_files "${source}/${pattern}") 
  install_("${src_files}" ${source} ${install_prefix})
endmacro()

macro(install_glob_recurse_if_newer_flat source pattern install_prefix)
  file(GLOB_RECURSE src_files "${source}/${pattern}") 
  install_flat_("${src_files}" ${source} ${install_prefix})
endmacro()

macro(install_glob_if_newer source pattern install_prefix)
  file(GLOB src_files "${source}/${pattern}") 
  install_("${src_files}" ${source} ${install_prefix})
endmacro()

## Use cmake_static_prefix, etc.
set(g_patterns "*.okl" "*.c" "*.hpp" "*.tpp" "*.h" "*.so" "*.dylib" "*.a" "hex*.dat")
set(g_suffix "okl" "include" "lib" "nodes" "data")

macro(install_batch_ source_dir sub_dir dest_dir patterns)
  foreach(pat ${patterns})
    install_glob_if_newer(${source_dir} ${pat} ${dest_dir})
    foreach(sub ${sub_dir})
      install_glob_recurse_if_newer(${source_dir}/${sub} ${pat} ${dest_dir}/${sub})
    endforeach()
  endforeach()
endmacro()

macro(install_batch_if_newer source_dir dest_dir)
    install_batch_(${source_dir} "${g_suffix}" ${dest_dir} "${g_patterns}")
endmacro()

macro(install_files_if_newer source_files dest_dir)
  foreach(source ${source_files})
    get_filename_component(fname ${source} NAME)
    execute_process(
      COMMAND bash -c "mkdir -p ${dest_dir}"
    )
    execute_process(
      COMMAND bash -c "[ ${source} -nt ${dest_dir}/${fname} ] && cp ${source} ${dest_dir}/${fname}"
    )
  endforeach()
endmacro()

###############################################################################
# Install                                                                     #
###############################################################################
## libParanumal okl, headers, libs and nodes directory
message("-- Installing libParanumal")
install_batch_if_newer(@LIBPDIR@ @CMAKE_INSTALL_PREFIX@/libparanumal)
install_batch_if_newer(@OGSDIR@ @CMAKE_INSTALL_PREFIX@/gatherScatter)
install_batch_if_newer(@PARALMONDDIR@ @CMAKE_INSTALL_PREFIX@/parAlmond)
install_batch_if_newer(@ELLIPTICDIR@ @CMAKE_INSTALL_PREFIX@/elliptic)

## nek5000
message("-- Installing nek5000")
install_glob_recurse_if_newer(@NEKDIR@/core * @CMAKE_INSTALL_PREFIX@/nek5000/core)

install_batch_if_newer(@NEKDIR@/3rd_party/gslib
  @CMAKE_INSTALL_PREFIX@/nek5000/3rd_party/gslib)
install_batch_if_newer(@NEKDIR@/3rd_party/gslib
  @CMAKE_INSTALL_PREFIX@/nek5000/3rd_party/gslib)
install_batch_if_newer(@NEKDIR@/3rd_party/blasLapack
  @CMAKE_INSTALL_PREFIX@/nek5000/3rd_party/blasLapack)
install_batch_if_newer(@NEKDIR@/3rd_party/parRSB
  @CMAKE_INSTALL_PREFIX@/nek5000/3rd_party/parRSB)
install_batch_if_newer(@NEKDIR@/3rd_party/parRSB
  @CMAKE_INSTALL_PREFIX@/nek5000/3rd_party/parRSB)
install_batch_if_newer(@NEKDIR@/3rd_party/parMETIS
  @CMAKE_INSTALL_PREFIX@/nek5000/3rd_party/parMETIS)

install_files_if_newer(@NEKDIR@/bin/nekconfig @CMAKE_INSTALL_PREFIX@/nek5000/bin)

## nekRS
message("-- Installing nekRS")
install_glob_if_newer(@CMAKE_SOURCE_DIR@/scripts * @CMAKE_INSTALL_PREFIX@/bin)
install_glob_recurse_if_newer(@CMAKE_SOURCE_DIR@/src "*.h" @CMAKE_INSTALL_PREFIX@/include)
install_glob_recurse_if_newer(@CMAKE_SOURCE_DIR@/src "*.hpp" @CMAKE_INSTALL_PREFIX@/include)
install_glob_recurse_if_newer(@CMAKE_SOURCE_DIR@/okl "*.okl" @CMAKE_INSTALL_PREFIX@/okl)
#install_glob_recurse_if_newer(@CMAKE_SOURCE_DIR@/src/plugins/include "*.hpp" @CMAKE_INSTALL_PREFIX@/include/plugins)

install_files_if_newer(@CMAKE_SOURCE_DIR@/src/udf/CMakeLists.txt @CMAKE_INSTALL_PREFIX@/udf)
install_files_if_newer("@CMAKE_SOURCE_DIR@/src/nekInterface/NEKINTF;@CMAKE_SOURCE_DIR@/src/nekInterface/nekInterface.f;@CMAKE_SOURCE_DIR@/src/nekInterface/Makefile" @CMAKE_INSTALL_PREFIX@/nekInterface)
install_glob_recurse_if_newer(@CMAKE_SOURCE_DIR@/examples * @CMAKE_INSTALL_PREFIX@/examples)
