install(PROGRAMS bp5dbg.py
  RENAME bp5dbg
  DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_scripts-runtime
)
install(
  FILES
    adios2/bp5dbg/__init__.py
    adios2/bp5dbg/utils.py
    adios2/bp5dbg/metadata.py
    adios2/bp5dbg/metametadata.py
    adios2/bp5dbg/idxtable.py
  DESTINATION ${CMAKE_INSTALL_PYTHONDIR}/adios2/bp5dbg COMPONENT adios2_scripts-runtime
)
