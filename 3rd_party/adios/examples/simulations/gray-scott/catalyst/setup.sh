# The ${PARAVIEW_BUILD_PATH}/lib directory should contain `libADIOSInSituPlugin.so`.
export PARAVIEW_BUILD_PATH=/home/adios/Software/paraview/build
# Set the following env variables.
export ADIOS2_PLUGIN_PATH=${PARAVIEW_BUILD_PATH}/lib
export CATALYST_IMPLEMENTATION_NAME=paraview
export CATALYST_IMPLEMENTATION_PATHS=${PARAVIEW_BUILD_PATH}/lib/catalyst
