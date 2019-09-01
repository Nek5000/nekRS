#if !defined(nekrs_libparanumal_hpp_)
#define nekrs_libparanumal_hpp_

#include "mesh.h"
#include "mesh3D.h"
#include "ins.h"

namespace libParanumal {
  // Data structures
  using ::mesh_t;
  using ::ins_t;
  using ::mesh3D;
  using ::setupAide;

  // Functions
  using ::readDfloatArray;
  using ::readIntArray;
  // mesh
  using ::meshParallelConnect;
  using ::meshHaloSetup;
  using ::meshConnectBoundary;
  using ::meshLoadReferenceNodesHex3D;
  using ::meshConnectFaceNodes3D;
  using ::meshParallelConnectNodes;
  using ::meshSurfaceGeometricFactorsHex3D;
  using ::meshGeometricFactorsHex3D;
  // ins Solver
  using ::insSetup;
  using ::insError;
  using ::insGradient;
  using ::insVelocityRhs;
  using ::insVelocitySolve;
  using ::insVelocityUpdate;
  using ::insPressureRhs;
  using ::insPressureSolve;
  using ::insPressureUpdate;
  using ::insSubCycle;
  using ::insStrongSubCycle;
  using ::insAdvection;
  using ::insRestartWrite;
  using ::insReport;
}

#endif
