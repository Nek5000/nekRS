#if !defined(remap_vector_hpp_)
#define remap_vector_hpp_

#include <occa.hpp>
#include "nrssys.hpp"
#include <vector>
#include <type_traits>

class mesh_t;

// Given a mesh with a given polynomial order,
// construct mappings from E to L and L to E vector formats
template<typename FPType>
class LVector_t
{
public:
  
  static_assert(std::is_floating_point<FPType>::value, "FPType must be a floating point type");
  
  // at the moment, only dfloat and pfloat are supported
  static_assert(std::is_same<FPType,dfloat>::value || std::is_same<FPType,pfloat>::value, "FPType must be either dfloat or pfloat");
  
  LVector_t(const std::vector<mesh_t*>& meshes, bool oallocate = true);

  LVector_t(const LVector_t&) = delete;
  LVector_t& operator=(const LVector_t&) = delete;
  
  occa::memory & optr();
  const occa::memory & optr() const;
  void optr(const std::vector<dlong> &fieldOffsets, occa::memory &o_L);
  void optr(occa::memory & o_L);

  void fieldOffsets(const std::vector<dlong> &fieldOffsets);
  const std::vector<dlong> &fieldOffsets() const;
  const std::vector<dlong> &fieldOffsetScans() const;

  const occa::memory & invDegree(int field) const;

  void copyToE(const dlong EFieldOffset, occa::memory & o_E) const;
  void copyFromE(const dlong EFieldOffset, const occa::memory & o_E);

  dlong nFields() const;
  dlong fieldOffset(int field) const;
  dlong fieldOffsetScan(int field) const;
  dlong Nlocal(int field) const;
  dlong size() const;

  const std::vector<mesh_t*> & meshes() const;

  static void registerKernels();

private:
  
  // convenience function to map FPType to string
  static std::string FPTypeString();
  
  static constexpr dlong UNMAPPED = -1;
  
  const std::vector<mesh_t*> meshes_;
  bool oallocate_;

  std::vector<dlong> fieldOffset_;
  std::vector<dlong> fieldOffsetScan_;
  std::vector<dlong> Nlocals_;
  dlong fieldOffsetSum_;

  occa::kernel EToLKernel_;
  occa::kernel LToEKernel_;

  occa::memory o_L_;

  std::map<mesh_t*, occa::memory> o_EToLs_;
  std::map<mesh_t*, occa::memory> o_EToLUniques_;
  std::map<mesh_t*, occa::memory> o_invDegrees_;

  std::map<mesh_t*, dlong> NfieldsMesh_; // Number of fields on each mesh
  std::map<mesh_t*, dlong> NlocalMesh_; // Nlocal on each mesh
  std::map<mesh_t*, occa::memory> fields_; // Map from field number inside mesh to actual field number
  
  // field information
  occa::memory o_fieldOffsetScan_;
};

#endif