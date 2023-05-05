#if !defined(nekrs_lpm_hpp_)
#define nekrs_lpm_hpp_

#include "nrssys.hpp"
#include <functional>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include "pointInterpolation.hpp"

class nrs_t;

// Lagrangian particle manager
class lpm_t {
public:
  using VerbosityLevel = pointInterpolation_t::VerbosityLevel;

  // User-defined RHS function
  // INPUTS:
  //   nrs: NekRS object
  //   lpm: particle manager
  //   time: current time
  //   o_y: particle degrees of freedom
  //        DO NOT ALTER!
  //   userdata: user-defined data
  // OUTPUTS:
  //   o_ydot: derivatives
  using rhsFunc_t = std::function<
      void(nrs_t *nrs, lpm_t *lpm, dfloat time, occa::memory o_y, void *userdata, occa::memory o_ydot)>;

  // User-defined ODE solver
  // Integrate from t0 to tf
  // INPUTS:
  //   nrs: NekRS object
  //   lpm: particle manager
  //   time: current time
  //   dt: time step(s)
  //   step: time step number
  //   o_y: particle degrees of freedom
  //   userdata: user-defined data
  //   o_ydot: scratch space
  using odeSolverFunc_t = std::function<void(nrs_t *nrs,
                                             lpm_t *lpm,
                                             dfloat t0,
                                             dfloat tf,
                                             int step,
                                             occa::memory o_y,
                                             void *userdata,
                                             occa::memory o_ydot)>;

  lpm_t(nrs_t *nrs, dfloat bb_tol_ = 0.01, dfloat newton_tol_ = 0.0);

  ~lpm_t() = default;

  // set AB integration order
  void abOrder(int order);

  // set RK integration order (1,2,3 or 4)
  void rkOrder(int order);

  // set ODE solver
  void setSolver(std::string solver);

  // static kernel registration function
  // This *MUST* be called during UDF_LoadKernels
  static void registerKernels(occa::properties &kernelInfo);

  // Register a degree of freedom (case insensitive)
  // By default, the particle coordinates are already registered
  // An optional bool argument may be used to flag a field as not needing to be output
  // during a lpm_t::writeFld(...) call.
  // By default, only the particle coordinates are automatically registered.
  // Pre:
  //   initialized() = false
  void registerDOF(std::string dofName, bool output = true);

  // Multi-component field version
  // On output, this field will be output to the VTU file as a vector quantity.
  void registerDOF(dlong Nfields, std::string dofName, bool output = true);

  // Properties associated with the particle
  // NOTE: these are _not_ degrees of freedom. Use registerDOF(...) for that.
  // As registerDOF(...), there's an optional bool output argument to
  // flag whether a field should be output during a lpm_t::writeFld(...) call.
  // Pre:
  //   initialized() = false
  void registerProp(std::string propName, bool output = true);

  // Multi-component field version
  // On output, this field will be output to the VTU file as a vector quantity.
  void registerProp(dlong Nfields, std::string propName, bool output = true);

  // Fields associated with the fluid mesh to be interpolated
  // to the particle locations.
  // On input, o_fld is a field on the fluid mesh to be interpolated.
  // NOTE: these are _not_ degrees of freedom. Use registerDOF(...) for that.
  // As registerDOF(...), there's an optional bool output argument to
  // flag whether a field should be output during a lpm_t::writeFld(...) call.
  // Pre:
  //   initialized() = false
  void registerInterpField(std::string interpFieldName, occa::memory o_fld, bool output = true);

  // Multi-component field version
  // Prefer using this version for multi-component fields, as the performance
  // will be better during interpolation.
  // On output, this field will be output to the VTU file as a vector quantity.
  void
  registerInterpField(std::string interpFieldName, dlong Nfields, occa::memory o_fld, bool output = true);

  // Get field index associated with a degree of freedom
  int dofId(std::string dofName) const;

  // Number of fields associated with a DOF
  int numDOFs(std::string dofName) const;

  // Get field index associated with a property
  int propId(std::string propName) const;

  // Number of fields associated with a property
  int numProps(std::string propName) const;

  // Get field index associated with an interpolated field
  int interpFieldId(std::string interpFieldName) const;

  // Number of fields associated with an interpolated field
  int numFieldsInterp(std::string interpFieldName) const;

  bool initialized() const { return initialized_; }

  // Required user RHS
  void setUserRHS(rhsFunc_t userRHS);

  // Optionally set user ODE solver
  void setUserODESolver(odeSolverFunc_t userODESolver) { userODESolver_ = userODESolver; };

  // Add optional userdata ptr to be passed to userRHS
  void addUserData(void *userdata);

  // Create nParticles particles with initial condition y0
  void initialize(int nParticles, dfloat t0, std::vector<dfloat> &y0);

  // Create nParticles particles with initial condition o_y0
  void initialize(int nParticles, dfloat t0, occa::memory o_y0);

  // Page-aligned offset >= nParticles
  // Required to access particle-specific fields
  // NOT valid until construct() is called
  int fieldOffset() const { return fieldOffset_; }

  // Compute page-aligned offset >= n
  static int fieldOffset(int n);

  // Number of particle degrees of freedom
  int nDOFs() const { return nDOFs_; }

  // Number of particle properties
  int nProps() const { return nProps_; }

  // Number of interpolated fields
  int nInterpFields() const { return nInterpFields_; }

  // Number of particles
  int numParticles() const { return nParticles_; }

  // Integrate to state tf
  // Pre:
  //   initialized() = true
  void integrate(dfloat tf);

  // Write particle data to file
  void writeFld();

  // Read particle data from file
  // Can be called in lieu of construct
  void restart(std::string restartFile);

  // Get particle degrees of freedom on device
  // Pre:
  //  initialized() = true
  occa::memory getDOF(int dofId);
  occa::memory getDOF(std::string dofName);

  // Get particle coordinates on host
  std::vector<dfloat> getDOFHost(std::string dofName);

  // Get particle property on device
  // Pre:
  //  initialized() = true
  occa::memory getProp(int propId);
  occa::memory getProp(std::string propName);

  // Get particle properties on host
  std::vector<dfloat> getPropHost(std::string propName);

  // Get interpolated field on device
  // Pre:
  //  initialized() = true
  occa::memory getInterpField(int interpFieldId);
  occa::memory getInterpField(std::string interpFieldName);

  // Get interpolated fields on host
  std::vector<dfloat> getInterpFieldHost(std::string interpFieldName);

  // Get the underlying pointInterpolation_t object
  pointInterpolation_t &interpolator() { return *interp; }

  // Interpolate all interpolated fields from the fluid mesh to the particle locations
  void interpolate();

  // Interpolate specific field name from the fluid mesh to the particle locations
  void interpolate(std::string interpFieldName);

  // Add new particles
  // The lagged ydot values are assumed to be zero in this case --
  // if accuracy matters, consider using an explicit RK method as the integrator,
  // or call the addParticles(...) overload that takes ydot as an argument.
  void addParticles(int nParticles, occa::memory o_yNewPart, occa::memory o_propNewPart);
  void
  addParticles(int nParticles, const std::vector<dfloat> &yNewPart, const std::vector<dfloat> &propNewPart);

  void addParticles(int nParticles,
                    occa::memory o_yNewPart,
                    occa::memory o_propNewPart,
                    occa::memory o_ydotNewPart);
  void addParticles(int nParticles,
                    const std::vector<dfloat> &yNewPart,
                    const std::vector<dfloat> &propNewPart,
                    const std::vector<dfloat> &ydotNewPart);

  // Delete particles that have left the domain
  void deleteParticles();

  // Moves particles to the processor able to evaluate them
  void migrate();

  // Number of particles across all MPI ranks
  long long int numGlobalParticles() const;

  // Number of particles that cannot be found inside the fluid mesh
  int numUnfoundParticles() const;

  // Number of particles that cannot be processed on the current rank
  int numNonLocalParticles() const;

  occa::memory o_prop;      // particle properties
  occa::memory o_interpFld; // interpolated field outputs
  occa::memory o_y;         // particle degrees of freedom
  occa::memory o_ydot;      // derivatives

  // set timer level
  void setTimerLevel(TimerLevel level);
  TimerLevel getTimerLevel() const;

  // set timer name
  // this is used to prefix the timer names
  void setTimerName(std::string name);

  // enables verbose warnings from findpts
  void setVerbosity(VerbosityLevel level) { verbosityLevel = level; };

  void printTimers();

  void resetTimers();

private:
  dlong nEXT, nBDF;

  static constexpr int bootstrapRKOrder = 4;
  
  // maximum number of entries valid for lpm_t::migration call
  static constexpr int maxEntriesPerParticleMigration = 50;

  enum class SolverType { AB, RK, INVALID };

  occa::memory o_ytmp; // scratch memory for RK integrators
  occa::memory o_k;    // k1, ... for RK integrators

  static SolverType stringToSolverType(std::string solverType);

  SolverType solverType = SolverType::AB;

  // Required to handle runtime -> compile time switch needed for sarray_transfer
  template <int N>
  void sendReceiveDataImpl(const std::vector<dfloat> &sendData,
                           const std::vector<dfloat> &r,
                           const std::vector<dlong> &proc,
                           const std::vector<dlong> &code,
                           const std::vector<dlong> &elem,
                           std::vector<dfloat> &receivedData,
                           std::vector<dfloat> &recvR,
                           std::vector<dlong> &recvCode,
                           std::vector<dlong> &recvElem,
                           int entriesPerParticle);

  // helper to return results after a migration call
  std::tuple<std::vector<dfloat>, std::vector<dfloat>, std::vector<dlong>, std::vector<dlong>>
  sendReceiveData(const std::vector<dfloat> &sendData,
                  const std::vector<dfloat> &r,
                  const std::vector<dlong> &proc,
                  const std::vector<dlong> &code,
                  const std::vector<dlong> &elem,
                  int entriesPerParticle);

  // helper function to handle allocations dependent on nParticles
  void handleAllocation(int offset);

  // helper function to call findpts to set up interpolation, given a new particle locations
  void find(occa::memory o_yNew);

  // helper function to extrapolate fluid state to a specified time
  void extrapolateFluidState(dfloat tEXT);

  // implements AB integrator
  void integrateAB();

  // implements RK integrator
  void integrateRK();
  void integrateRK1();
  void integrateRK2();
  void integrateRK3();
  void integrateRK4();

  // generate set of all output DOFs, sans {x,y,z}
  std::set<std::string> nonCoordinateOutputDOFs() const;

  void abCoeff(dfloat *dt, int tstep);

  dfloat time = 0.0;
  int tstep = 0;

  std::string timerName = "lpm_t::";
  TimerLevel timerLevel = TimerLevel::Basic;
  VerbosityLevel verbosityLevel = VerbosityLevel::Basic;
  nrs_t *nrs = nullptr;
  int solverOrder;
  dfloat newton_tol;
  dfloat bb_tol;
  std::unique_ptr<pointInterpolation_t> interp;

  int nParticles_ = 0;
  int nDOFs_ = 0;
  int nProps_ = 0;
  int nInterpFields_ = 0;
  int fieldOffset_ = 0; // page-aligned offset >= nParticles
  bool initialized_ = false;
  inline static bool kernelsRegistered_ = false;

  std::vector<dfloat> coeffAB;
  occa::memory o_coeffAB; // AB coefficients

  std::vector<dfloat> coeffRK;
  occa::memory o_coeffRK; // RK coefficients

  std::vector<dfloat> coeffEXT;
  occa::memory o_coeffEXT;   // EXT coefficients
  std::vector<dfloat> dtEXT; // for extrapolating velocity states

  std::vector<dfloat> dt; // previous time steps

  rhsFunc_t userRHS_ = nullptr;
  odeSolverFunc_t userODESolver_ = nullptr;

  enum class FieldType { DOF, PROP, INTERP_FIELD };

  // Map between dof/prop/interpField names and which type of field they are
  // This is especially useful when reading in restart files
  std::map<std::string, FieldType> fieldType;

  // DOFs
  std::map<std::string, int> dofIds;
  std::map<std::string, int> dofCounts;
  std::map<std::string, bool> outputDofs;

  // Properties
  std::map<std::string, int> propIds;
  std::map<std::string, int> propCounts;
  std::map<std::string, bool> outputProps;

  // Interpolated fields
  std::map<std::string, int> interpFieldIds;
  std::map<std::string, int> interpFieldCounts;
  std::map<std::string, bool> outputInterpFields;
  std::map<std::string, occa::memory> interpFieldInputs;

  // History of interpolated fields
  std::map<std::string, occa::memory> laggedInterpFields;

  // Hold extrapolated state during a particle integration
  std::map<std::string, occa::memory> extrapolatedInterpFields;

  // map from current particle id to new particle id when adding/deleting particles
  occa::memory o_remainingMap;

  // map new particles to particle id when adding particles
  occa::memory o_insertMap;

  // map from current particle id to new particle id when migrating particles
  occa::memory o_currentRankMap;

  // map from current particle id to sending rank when migrating particles
  occa::memory o_sendRankMap;

  // map new particles to particle id when migrating particles
  occa::memory o_migrateMap;

  // map from received particle to particle id when migrating particles
  occa::memory o_recvRankMap;

  void *userdata_ = nullptr;
  occa::kernel nStagesSumManyKernel;
  occa::kernel remapParticlesKernel;
};

#endif
