#include "LVector.hpp"
#include "mesh.h"
#include "platform.hpp"
#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace
{

struct LVectorMapping_t {
  dlong Nlocal;
  occa::memory o_EToL;
  occa::memory o_EToLUnique;
  occa::memory o_invDegree;
  occa::memory o_invDegreePfloat;
};

// Handle allocations
class LVectorMappingManager_t
{
public:
  static LVectorMappingManager_t &instance()
  {
    static LVectorMappingManager_t instance_;
    return instance_;
  }

  LVectorMappingManager_t(const LVectorMappingManager_t &) = delete;
  LVectorMappingManager_t &operator=(const LVectorMappingManager_t &) = delete;

  auto get(mesh_t *mesh)
  {
    if (mappings_.count(mesh) == 0) {
      mappings_[mesh] = setup(mesh);
    }

    return mappings_[mesh];
  }

private:
  LVectorMappingManager_t() = default;

  LVectorMapping_t setup(mesh_t *mesh)
  {
    auto o_Lids = platform->deviceMemoryPool.reserve<dlong>(mesh->Nlocal);
    std::vector<dlong> Eids(mesh->Nlocal);
    std::iota(Eids.begin(), Eids.end(), 0);
    o_Lids.copyFrom(Eids.data(), mesh->Nlocal);

    {
      const auto saveNhaloGather = mesh->ogs->NhaloGather;
      mesh->ogs->NhaloGather = 0;
      ogsGatherScatter(o_Lids, ogsInt, ogsMin, mesh->ogs);
      mesh->ogs->NhaloGather = saveNhaloGather;
    }

    std::vector<dlong> Lids(mesh->Nlocal);
    o_Lids.copyTo(Lids.data(), mesh->Nlocal);

    std::set<dlong> uniqueIds;
    for (auto &&id : Lids) {
      uniqueIds.insert(id);
    }

    const int NL = uniqueIds.size();

    std::vector<dlong> EToLUnique(mesh->Nlocal, -1);
    unsigned ctr = 0;
    for (auto &&uniqueEid : uniqueIds) {
      EToLUnique[uniqueEid] = ctr;
      ctr++;
    }
    auto o_EToLUnique = platform->device.malloc<dlong>(mesh->Nlocal, EToLUnique.data());

    // setup non-unique version mapping an E-vector point to its corresponding L-vector point
    std::vector<dlong> EToL(mesh->Nlocal, -1);
    for (int n = 0; n < mesh->Nlocal; ++n) {
      const auto lid = Lids[n];
      EToL[n] = EToLUnique[lid];
    }
    auto o_EToL = platform->device.malloc<dlong>(mesh->Nlocal, EToL.data());

    // construct L-vector version of inv degree, based on duplicated points in L-vector
    auto [o_invDegreeL, o_invDegreeLPfloat] = [&]() {
      std::vector<dfloat> degree(mesh->Nlocal, 0.0);

      for (int n = 0; n < mesh->Nlocal; ++n) {
        const auto lid = EToLUnique[n];
        if (lid > -1) {
          degree[n] = 1.0;
        }
      }

      ogsGatherScatter(degree.data(), dfloatString, ogsAdd, mesh->ogs);

      std::vector<dfloat> invDegreeL(NL, 1.0);

      for (int n = 0; n < mesh->Nlocal; ++n) {
        const auto lid = EToLUnique[n];
        if (lid > -1) {
          invDegreeL[lid] = 1.0 / degree[n];
        }
      }

      // sanity check:
      // entries in invDegreeL are on (0,1]
      bool allPositive = std::all_of(invDegreeL.begin(), invDegreeL.end(), [](auto &&val) {
        return val > 0.0 && val <= 1.0;
      });
      int err = allPositive ? 0 : 1;
      MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
      nekrsCheck(err,
                 platform->comm.mpiComm,
                 EXIT_FAILURE,
                 "%s\n",
                 "Encountered invDegreeL value outside of (0,1]");

      // on a single processor, T-vector and L-vector are the same --> invDegreeL is unity everywhere
      if (platform->comm.mpiCommSize == 1) {
        const auto tol = 1e4 * std::numeric_limits<dfloat>::epsilon();
        bool allUnity = std::all_of(invDegreeL.begin(), invDegreeL.end(), [tol](auto &&val) {
          return std::abs(val - 1.0) < tol;
        });
        int err = allUnity ? 0 : 1;
        MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
        nekrsCheck(err,
                   platform->comm.mpiComm,
                   EXIT_FAILURE,
                   "%s\n",
                   "Encountered non-unity invDegreeL when P=1");
      }

      // construct pfloat version
      std::vector<pfloat> invDegreeLPfloat(NL, 1.0);
      for (int n = 0; n < NL; ++n) {
        invDegreeLPfloat[n] = static_cast<pfloat>(invDegreeL[n]);
      }

      auto o_invDegreeL = platform->device.malloc<dfloat>(NL, invDegreeL.data());
      auto o_invDegreeLPfloat = platform->device.malloc<pfloat>(NL, invDegreeLPfloat.data());
      return std::make_tuple(o_invDegreeL, o_invDegreeLPfloat);
    }();

    // a few sanity checks:
    // EToL has only non-negative values
    bool allNonNegative = std::all_of(EToL.begin(), EToL.end(), [](auto &&val) { return val >= 0; });
    int err = allNonNegative ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    nekrsCheck(err,
               platform->comm.mpiComm,
               EXIT_FAILURE,
               "%s\n",
               "Encountered negative value in EToL mapping");

    // range of EToL is [0, NL)
    auto minmax = std::minmax_element(EToL.begin(), EToL.end());
    err = (*minmax.first == 0 && *minmax.second == NL - 1) ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "EToL mapping is not in range [0, NL)");

    // every value in [0,NL) is covered in the range of EToL
    std::set<dlong> uniqueOutputs;
    for (auto &&val : EToL) {
      uniqueOutputs.insert(val);
    }
    err = (uniqueOutputs.size() == NL) ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "EToL mapping is not surjective");

    o_Lids.free();

    return LVectorMapping_t{NL, o_EToL, o_EToLUnique, o_invDegreeL, o_invDegreeLPfloat};
  }

  std::map<mesh_t *, LVectorMapping_t> mappings_;
};

} // namespace

template <typename FPType>
LVector_t<FPType>::LVector_t(const std::vector<mesh_t *> &meshes, bool oallocate)
    : meshes_(meshes), oallocate_(oallocate), fieldOffset_(meshes.size(), 0),
      fieldOffsetScan_(meshes.size(), 0), Nlocals_(meshes.size(), 0)
{
  if (meshes.size() == 0) {
    return; // empty LVector
  }

  std::map<mesh_t *, std::vector<dlong>> meshToFields;
  for (auto &&mesh : meshes) {
    NfieldsMesh_[mesh] = 0;
    meshToFields[mesh] = std::vector<dlong>();
  }

  fieldOffsetSum_ = 0;
  for (int fld = 0; fld < meshes.size(); ++fld) {
    auto mesh = meshes[fld];
    auto mapping = LVectorMappingManager_t::instance().get(mesh);
    o_EToLs_[mesh] = mapping.o_EToL;
    o_EToLUniques_[mesh] = mapping.o_EToLUnique;

    if (std::is_same<FPType, dfloat>::value) {
      o_invDegrees_[mesh] = mapping.o_invDegree;
    } else {
      o_invDegrees_[mesh] = mapping.o_invDegreePfloat;
    }

    Nlocals_[fld] = mapping.Nlocal;
    fieldOffset_[fld] = alignStride<FPType>(mapping.Nlocal);
    fieldOffsetSum_ += fieldOffset_[fld];
    NfieldsMesh_[mesh]++;
    NlocalMesh_[mesh] = mapping.Nlocal;
    meshToFields[mesh].push_back(fld);
  }

  for (auto &&mesh : meshes) {
    auto fields = meshToFields.at(mesh);
    fields_[mesh] = platform->device.malloc<dlong>(fields.size(), fields.data());
  }

  fieldOffsetScan_[0] = 0;
  for (int fld = 1; fld < meshes.size(); ++fld) {
    fieldOffsetScan_[fld] = fieldOffsetScan_[fld - 1] + fieldOffset_[fld - 1];
  }

  o_fieldOffsetScan_ = platform->device.malloc<dlong>(meshes.size(), fieldOffsetScan_.data());
  if (oallocate) {
    o_L_ = platform->device.malloc<FPType>(fieldOffsetSum_);
  }

  auto suffix = LVector_t<FPType>::FPTypeString();
  this->EToLKernel_ = platform->kernelRequests.load("EToL" + suffix);
  this->LToEKernel_ = platform->kernelRequests.load("LToE" + suffix);
}

template <typename FPType> occa::memory &LVector_t<FPType>::optr()
{
  return this->o_L_;
}

template <typename FPType> const occa::memory &LVector_t<FPType>::optr() const
{
  return this->o_L_;
}

template <typename FPType>
void LVector_t<FPType>::optr(const std::vector<dlong> &fieldOffsets, occa::memory &o_L)
{
  this->fieldOffsets(fieldOffsets);
  this->optr(o_L);
}

template <typename FPType> void LVector_t<FPType>::optr(occa::memory &o_L)
{
  nekrsCheck(o_L.byte_size() <= this->size(),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "LVector_t::optr o_L.byte_size() = %llu, while expecting at least %d entries!\n",
             o_L.byte_size(),
             this->size());

  if (this->o_L_.isInitialized() && oallocate_) {
    this->o_L_.free();
  }

  this->o_L_ = o_L;
}

template <typename FPType> void LVector_t<FPType>::fieldOffsets(const std::vector<dlong> &fieldOffsets)
{
  if (this->nFields() == 0) {
    return;
  }

  // offsets size _must_ match the number of fields
  nekrsCheck(fieldOffsets.size() != this->nFields(),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "LVector_t::offsets offsets.size() = %ld, while expecting %d entries!\n",
             fieldOffsets.size(),
             this->nFields());

  // check that each fieldOffset is at least as large as the corresponding Nlocal
  std::ostringstream errLogger;
  for (int fld = 0; fld < this->nFields(); ++fld) {
    if (fieldOffsets[fld] < this->Nlocal(fld)) {
      errLogger << "LVector_t::offsets offsets[" << fld << "] = " << fieldOffsets[fld]
                << ", while expecting at least " << this->Nlocal(fld) << "!\n";
    }
  }

  const auto errString = errLogger.str();
  nekrsCheck(errString.size() != 0, platform->comm.mpiComm, EXIT_FAILURE, "%s", errString.c_str());

  std::copy(fieldOffsets.begin(), fieldOffsets.end(), this->fieldOffset_.begin());
  this->fieldOffsetScan_[0] = 0;
  this->fieldOffsetSum_ = this->fieldOffset_[0];
  for (int fld = 1; fld < this->nFields(); ++fld) {
    this->fieldOffsetScan_[fld] = this->fieldOffsetScan_[fld - 1] + this->fieldOffset_[fld - 1];
    this->fieldOffsetSum_ += this->fieldOffset_[fld];
  }
  this->o_fieldOffsetScan_.copyFrom(this->fieldOffsetScan_.data());
}

template <typename FPType> const std::vector<dlong> &LVector_t<FPType>::fieldOffsets() const
{
  return this->fieldOffset_;
}

template <typename FPType> const std::vector<dlong> &LVector_t<FPType>::fieldOffsetScans() const
{
  return this->fieldOffsetScan_;
}

template <typename FPType> const occa::memory &LVector_t<FPType>::invDegree(int field) const
{
  return this->o_invDegrees_.at(meshes_.at(field));
}

template <typename FPType> void LVector_t<FPType>::copyToE(const dlong EFieldOffset, occa::memory &o_E) const
{
  // loop over unique mesh objects to minimize number of kernel calls
  for (auto [uniqueMesh, _] : o_EToLs_) {
    LToEKernel_(uniqueMesh->Nlocal,
                NfieldsMesh_.at(uniqueMesh),
                EFieldOffset,
                fields_.at(uniqueMesh),
                o_fieldOffsetScan_,
                o_EToLs_.at(uniqueMesh),
                o_L_,
                o_E);
  }
}

template <typename FPType>
void LVector_t<FPType>::copyFromE(const dlong EFieldOffset, const occa::memory &o_E)
{
  // loop over unique mesh objects to minimize number of kernel calls
  for (auto [uniqueMesh, _] : o_EToLs_) {
    EToLKernel_(uniqueMesh->Nlocal,
                NfieldsMesh_.at(uniqueMesh),
                EFieldOffset,
                fields_.at(uniqueMesh),
                o_fieldOffsetScan_,
                o_EToLUniques_.at(uniqueMesh),
                o_E,
                o_L_);
  }
}

template <typename FPType> dlong LVector_t<FPType>::nFields() const
{
  return this->meshes_.size();
}

template <typename FPType> dlong LVector_t<FPType>::fieldOffset(int field) const
{
  return this->fieldOffset_.at(field);
}

template <typename FPType> dlong LVector_t<FPType>::fieldOffsetScan(int field) const
{
  return this->fieldOffsetScan_.at(field);
}

template <typename FPType> dlong LVector_t<FPType>::Nlocal(int field) const
{
  return this->Nlocals_.at(field);
}

template <typename FPType> dlong LVector_t<FPType>::size() const
{
  return this->fieldOffsetSum_;
}

template <typename FPType> const std::vector<mesh_t *> &LVector_t<FPType>::meshes() const
{
  return this->meshes_;
}

template <typename FPType> void LVector_t<FPType>::registerKernels()
{
  const std::string oklpath = getenv("NEKRS_KERNEL_DIR") + std::string("/core/");

  auto kernelInfo = platform->kernelInfo;
  kernelInfo["defines/"
             "p_UNMAPPED"] = UNMAPPED;

  if (std::is_same<FPType, dfloat>::value) {
    kernelInfo["defines/"
               "pfloat"] = "dfloat";
  } else {
    kernelInfo["defines/"
               "pfloat"] = "pfloat";
  }

  auto suffix = LVector_t<FPType>::FPTypeString();

  {
    auto kernelName = "EToL";
    auto fileName = oklpath + kernelName + ".okl";
    platform->kernelRequests.add(kernelName + suffix, fileName, kernelInfo);
  }

  {
    auto kernelName = "LToE";
    auto fileName = oklpath + kernelName + ".okl";
    platform->kernelRequests.add(kernelName + suffix, fileName, kernelInfo);
  }
}

template <typename FPType> std::string LVector_t<FPType>::FPTypeString()
{
  if (std::is_same<FPType, dfloat>::value) {
    return "dfloat";
  } else {
    return "pfloat";
  }
}

// instantiate dfloat, pfloat versions of LVector_t
template class LVector_t<dfloat>;
#ifndef NEKRS_USE_DFLOAT_FLOAT
template class LVector_t<pfloat>;
#endif
