install(PROGRAMS adios2_campaign_manager.py
  RENAME adios2_campaign_manager
  DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_scripts-runtime
)
#install(
#  FILES
#    adios2/bp4dbg/__init__.py
#    adios2/bp4dbg/data.py
#    adios2/bp4dbg/utils.py
#    adios2/bp4dbg/metadata.py
#    adios2/bp4dbg/idxtable.py
#  DESTINATION ${CMAKE_INSTALL_PYTHONDIR}/adios2/bp4dbg COMPONENT adios2_scripts-runtime
#)
