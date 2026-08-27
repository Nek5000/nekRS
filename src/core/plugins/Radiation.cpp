#include "nrs.hpp"
#include "platform.hpp"
#include "Radiation.hpp"
#include "RadiationBVH.hpp"
#include "mesh3D.h"
#include "sha1.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>

namespace
{

nrs_t *nrs;

bool buildKernelCalled = false;
bool setupCalled = false;

occa::kernel samplePairsKernel;

// Runtime state retained across Radiation::step() calls, populated once by
// Radiation::setup() and reused every RADIATION updateFrequency steps to
// couple the (fixed) view-factor matrix to the (time-varying) temperature
// field via a gray-diffuse radiosity solve.
int P_state = 0;
std::vector<dfloat> F_state;             // [P][P], rank 0 only
std::vector<dfloat> patchEmissivity_state; // [P], rank 0 only
std::vector<dfloat> J_state;              // [P] radiosity, rank 0 only, warm-started
dfloat stefanBoltzmann_state = (dfloat)5.670374419e-8;
int updateFrequency_state = 10;
dfloat radiosityTolerance_state = (dfloat)1e-6;
int radiosityMaxIters_state = 100;

// Boundary-ID groups (sorted ascending, matching the CSV output at setup
// time) and each patch's group index -- kept so Radiation::step() can print
// a per-group average flux each update, the same grouping used for the
// setup-time <case>_radiation_viewfactors_groups.csv.
std::vector<int> groups_state;
std::vector<int> patchGroupIdx_state; // [P]

// Which of this rank's local radiating patches (global patch index p, and
// the Nfp volume-node indices for that patch face) need a flux scattered
// back into bc->o_usrwrk each update.
struct LocalRadiatingPatch {
  int p;
  std::vector<dlong> idxVol;
};
std::vector<LocalRadiatingPatch> localRadiatingPatches_state;

std::vector<dfloat> parseDoubleList(const std::string &raw)
{
  std::vector<dfloat> vals;
  for (const auto &tok : serializeString(raw, ',')) {
    if (!tok.empty()) {
      vals.push_back((dfloat)std::stod(tok));
    }
  }
  return vals;
}

struct RadiationPatch {
  std::vector<dfloat> coords; // 3*Nfp, laid out coords[3*n+0/1/2]
  std::vector<dlong> idxVol;  // Nfp volume-node indices (mesh->vmapM), local to this rank
  int boundaryID = 0;
  dfloat refNx = 0, refNy = 0, refNz = 0;
};

std::vector<int> parseBoundaryIDList(const std::string &raw)
{
  std::vector<int> ids;
  for (const auto &tok : serializeString(raw, ',')) {
    if (!tok.empty()) {
      ids.push_back(std::stoi(tok));
    }
  }
  return ids;
}

// Gathers this rank's local boundary faces whose EToB matches one of the
// configured radiating/obstruction boundary IDs. isRadiatingFlag/
// isObstructionFlag are filled in parallel with the returned patch list.
std::vector<RadiationPatch> gatherLocalPatches(mesh_t *mesh,
                                               const std::set<int> &participatingIDs,
                                               const std::set<int> &radiatingIDs,
                                               const std::set<int> &obstructionIDs,
                                               std::vector<int> &isRadiatingFlag,
                                               std::vector<int> &isObstructionFlag)
{
  std::vector<RadiationPatch> patches;

  auto [x, y, z] = mesh->xyzHost();

  std::vector<dfloat> sgeoNode(mesh->Nsgeo);

  for (dlong e = 0; e < mesh->Nelements; ++e) {
    for (int f = 0; f < mesh->Nfaces; ++f) {
      const int bID = mesh->EToB[f + e * mesh->Nfaces];
      if (bID <= 0 || participatingIDs.find(bID) == participatingIDs.end()) {
        continue;
      }

      RadiationPatch patch;
      patch.boundaryID = bID;
      patch.coords.resize(3 * mesh->Nfp);
      patch.idxVol.resize(mesh->Nfp);

      for (int n = 0; n < mesh->Nfp; ++n) {
        const dlong idM = mesh->vmapM[e * mesh->Nfaces * mesh->Nfp + f * mesh->Nfp + n];
        patch.coords[3 * n + 0] = x[idM];
        patch.coords[3 * n + 1] = y[idM];
        patch.coords[3 * n + 2] = z[idM];
        patch.idxVol[n] = idM;
      }

      const dlong sid = e * mesh->Nfaces * mesh->Nfp + f * mesh->Nfp + 0;
      mesh->o_sgeo.copyTo(sgeoNode.data(), mesh->Nsgeo, sid * mesh->Nsgeo);
      // sgeo's normal points OUT of the fluid domain (CFD convention). A
      // radiating surface's normal must point INTO the enclosure it bounds
      // (toward whatever it can see), so negate it here -- empirically
      // confirmed via the radiationPlates case: with the un-negated normal,
      // cosI/cosJ were always negative for any cross-plate pair (both plates'
      // normals pointing away from each other), so every sample was rejected
      // and F_ij came out as exactly 0 for all off-diagonal pairs.
      patch.refNx = -sgeoNode[NXID];
      patch.refNy = -sgeoNode[NYID];
      patch.refNz = -sgeoNode[NZID];

      patches.push_back(patch);
      isRadiatingFlag.push_back(radiatingIDs.count(bID) ? 1 : 0);
      isObstructionFlag.push_back(obstructionIDs.count(bID) ? 1 : 0);
    }
  }

  return patches;
}

struct GlobalPatches {
  std::vector<dfloat> coords;    // [K][3*Nfp]
  std::vector<int> boundaryID;   // [K]
  std::vector<dfloat> refNormal; // [K][3]
  std::vector<int> isRadiating;  // [K]
  std::vector<int> isObstruction; // [K]
  int K = 0;
};

// Concatenates every rank's local patch list, in rank order, so every rank
// ends up with an identical global list and identical global indices for
// each physical face (used as the cross-rank "face key" for self-exclusion).
GlobalPatches allgatherPatches(const std::vector<RadiationPatch> &local,
                               const std::vector<int> &isRadiatingFlag,
                               const std::vector<int> &isObstructionFlag,
                               int Nfp,
                               MPI_Comm comm,
                               int &myDispl)
{
  const int nLocal = static_cast<int>(local.size());

  std::vector<dfloat> localCoords((size_t)nLocal * 3 * Nfp);
  std::vector<int> localBID(nLocal);
  std::vector<dfloat> localRefN((size_t)nLocal * 3);
  for (int i = 0; i < nLocal; ++i) {
    std::copy(local[i].coords.begin(), local[i].coords.end(), localCoords.begin() + (size_t)i * 3 * Nfp);
    localBID[i] = local[i].boundaryID;
    localRefN[3 * i + 0] = local[i].refNx;
    localRefN[3 * i + 1] = local[i].refNy;
    localRefN[3 * i + 2] = local[i].refNz;
  }

  int rank, nRanks;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);

  std::vector<int> counts(nRanks);
  MPI_Allgather(&nLocal, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

  std::vector<int> displs(nRanks);
  int total = 0;
  for (int r = 0; r < nRanks; ++r) {
    displs[r] = total;
    total += counts[r];
  }
  myDispl = displs[rank];

  GlobalPatches g;
  g.K = total;
  g.coords.resize((size_t)total * 3 * Nfp);
  g.boundaryID.resize(total);
  g.refNormal.resize((size_t)total * 3);
  g.isRadiating.resize(total);
  g.isObstruction.resize(total);

  std::vector<int> coordCounts(nRanks), coordDispls(nRanks);
  std::vector<int> vec3Counts(nRanks), vec3Displs(nRanks);
  for (int r = 0; r < nRanks; ++r) {
    coordCounts[r] = counts[r] * 3 * Nfp;
    coordDispls[r] = displs[r] * 3 * Nfp;
    vec3Counts[r] = counts[r] * 3;
    vec3Displs[r] = displs[r] * 3;
  }

  MPI_Allgatherv(localCoords.data(),
                nLocal * 3 * Nfp,
                MPI_DFLOAT,
                g.coords.data(),
                coordCounts.data(),
                coordDispls.data(),
                MPI_DFLOAT,
                comm);

  MPI_Allgatherv(localBID.data(), nLocal, MPI_INT, g.boundaryID.data(), counts.data(), displs.data(), MPI_INT, comm);

  MPI_Allgatherv(localRefN.data(),
                nLocal * 3,
                MPI_DFLOAT,
                g.refNormal.data(),
                vec3Counts.data(),
                vec3Displs.data(),
                MPI_DFLOAT,
                comm);

  MPI_Allgatherv(isRadiatingFlag.data(),
                nLocal,
                MPI_INT,
                g.isRadiating.data(),
                counts.data(),
                displs.data(),
                MPI_INT,
                comm);

  MPI_Allgatherv(isObstructionFlag.data(),
                nLocal,
                MPI_INT,
                g.isObstruction.data(),
                counts.data(),
                displs.data(),
                MPI_INT,
                comm);

  return g;
}

} // namespace

void Radiation::buildKernel(occa::properties _kernelInfo)
{
  occa::properties kernelInfo = _kernelInfo;

  auto build = [&kernelInfo](const std::string &kernelName) {
    const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/core/plugins/");
    const auto fileName = path + "Radiation.okl";
    const auto reqName = "Radiation::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = true;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  samplePairsKernel = build("radiationSamplePairs");
}

void Radiation::setup()
{
  static bool isInitialized = false;
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  if (!platform->options.compareArgs("RADIATION", "TRUE")) {
    return;
  }

  nekrsCheck(!buildKernelCalled,
            MPI_COMM_SELF,
            EXIT_FAILURE,
            "%s\n",
            "Radiation::setup called prior to Radiation::buildKernel!");

  nrs = dynamic_cast<nrs_t *>(platform->app);
  // meshV (not fluid->mesh) so this module works whether or not the flow
  // solver is enabled (FLUID=FALSE leaves nrs->fluid null) -- confirmed via
  // an empirical probe case that meshV is populated unconditionally.
  auto mesh = nrs->meshV;

  const int rank = platform->comm.mpiRank();
  MPI_Comm comm = platform->comm.mpiComm();

  std::string radiatingStr, obstructionStr;
  platform->options.getArgs("RADIATION RADIATING BOUNDARY IDS", radiatingStr);
  platform->options.getArgs("RADIATION OBSTRUCTION BOUNDARY IDS", obstructionStr);

  const auto radiatingList = parseBoundaryIDList(radiatingStr);
  const auto obstructionList = parseBoundaryIDList(obstructionStr);

  const std::set<int> radiatingIDs(radiatingList.begin(), radiatingList.end());
  const std::set<int> obstructionIDs(obstructionList.begin(), obstructionList.end());

  nekrsCheck(radiatingIDs.empty(), comm, EXIT_FAILURE, "%s\n", "[RADIATION] radiatingBoundaryIDs is empty!");

  std::set<int> participatingIDs = radiatingIDs;
  participatingIDs.insert(obstructionIDs.begin(), obstructionIDs.end());

  const bool hasObstruction = !obstructionIDs.empty();

  std::vector<int> isRadiatingFlag, isObstructionFlag;
  auto localPatches =
      gatherLocalPatches(mesh, participatingIDs, radiatingIDs, obstructionIDs, isRadiatingFlag, isObstructionFlag);

  int myDispl = 0;
  const auto global = allgatherPatches(localPatches, isRadiatingFlag, isObstructionFlag, mesh->Nfp, comm, myDispl);

  nekrsCheck(global.K == 0,
            comm,
            EXIT_FAILURE,
            "%s\n",
            "[RADIATION] no boundary faces matched the configured boundary IDs!");

  std::vector<int> radiatingIndices, obstructionIndices;
  for (int k = 0; k < global.K; ++k) {
    if (global.isRadiating[k]) {
      radiatingIndices.push_back(k);
    }
    if (global.isObstruction[k]) {
      obstructionIndices.push_back(k);
    }
  }

  const int P = static_cast<int>(radiatingIndices.size());
  nekrsCheck(P == 0, comm, EXIT_FAILURE, "%s\n", "[RADIATION] no radiating patches found!");

  // Output storage and the pair-sampling kernel both scale as O(P^2); this
  // module targets a modest "surface enclosure" subset of the mesh, not the
  // full CFD boundary. P=8000 -> ~32M packed pairs -> ~256MB for the raw
  // double matrix alone, already a reasonable practical ceiling for phase 1.
  constexpr int maxPatches = 8000;
  nekrsCheck(P > maxPatches,
            comm,
            EXIT_FAILURE,
            "[RADIATION] %d radiating patches exceeds the phase-1 limit of %d (O(P^2) memory/compute); "
            "narrow radiatingBoundaryIDs to a smaller enclosure.\n",
            P,
            maxPatches);

  if (rank == 0) {
    printf("Radiation: %d radiating patches, %d obstruction patches (global)\n",
          P,
          (int)obstructionIndices.size());
  }

  P_state = P;

  // Map each local patch (this rank's slice of the global, rank-order
  // concatenated patch list) back to its radiating-patch index p, so the
  // per-step flux update knows which of bc->o_usrwrk's nodes belong to which
  // patch. Radiating patches are not necessarily contiguous in p within a
  // rank's local range (an obstruction-only patch can fall in between), so
  // this is a lookup, not an offset.
  std::vector<int> globalKToP(global.K, -1);
  for (int p = 0; p < P; ++p) {
    globalKToP[radiatingIndices[p]] = p;
  }
  localRadiatingPatches_state.clear();
  for (int i = 0; i < static_cast<int>(localPatches.size()); ++i) {
    const int k = myDispl + i;
    const int p = globalKToP[k];
    if (p >= 0) {
      localRadiatingPatches_state.push_back({p, localPatches[i].idxVol});
    }
  }

  // Boundary-ID groups (same grouping used for the CSV output below), needed
  // on every rank -- not just rank 0 -- so the emissivity list can be parsed
  // and length-validated identically everywhere.
  std::set<int> groupIDSet;
  for (int p = 0; p < P; ++p) {
    groupIDSet.insert(global.boundaryID[radiatingIndices[p]]);
  }
  const std::vector<int> groups(groupIDSet.begin(), groupIDSet.end());
  const int nGroups = static_cast<int>(groups.size());

  std::vector<int> patchGroupIdx(P);
  for (int p = 0; p < P; ++p) {
    const int bID = global.boundaryID[radiatingIndices[p]];
    const int gi = static_cast<int>(std::lower_bound(groups.begin(), groups.end(), bID) - groups.begin());
    patchGroupIdx[p] = gi;
  }
  groups_state = groups;
  patchGroupIdx_state = patchGroupIdx;

  std::string emissivityStr;
  std::vector<dfloat> emissivity(nGroups, (dfloat)1.0);
  if (platform->options.getArgs("RADIATION EMISSIVITY", emissivityStr) && !emissivityStr.empty()) {
    emissivity = parseDoubleList(emissivityStr);
    nekrsCheck(static_cast<int>(emissivity.size()) != nGroups,
              comm,
              EXIT_FAILURE,
              "[RADIATION] emissivity has %d entries but there are %d radiating boundary-ID groups!\n",
              static_cast<int>(emissivity.size()),
              nGroups);
  }

  platform->options.getArgs("RADIATION STEFAN BOLTZMANN", stefanBoltzmann_state);
  platform->options.getArgs("RADIATION UPDATE FREQUENCY", updateFrequency_state);
  platform->options.getArgs("RADIATION RADIOSITY TOLERANCE", radiosityTolerance_state);
  platform->options.getArgs("RADIATION RADIOSITY MAX ITERS", radiosityMaxIters_state);

  if (rank == 0) {
    patchEmissivity_state.resize(P);
    for (int p = 0; p < P; ++p) {
      patchEmissivity_state[p] = emissivity[patchGroupIdx[p]];
    }
  }

  // Radiative flux BC storage: this module owns platform->app->bc->o_usrwrk
  // for phase 1 (single-purpose radiating cases). A case combining this with
  // another o_usrwrk-based BC would need a follow-up to share offsets.
  nekrsCheck(platform->app->bc->o_usrwrk.isInitialized() &&
                platform->app->bc->o_usrwrk.size() != (size_t)mesh->Nlocal,
            comm,
            EXIT_FAILURE,
            "%s\n",
            "[RADIATION] platform->app->bc->o_usrwrk is already sized for something else; "
            "the Radiation module currently requires exclusive use of it.");
  if (!platform->app->bc->o_usrwrk.isInitialized() || platform->app->bc->o_usrwrk.size() == 0) {
    platform->app->bc->o_usrwrk.resize(mesh->Nlocal);
    std::vector<dfloat> zeroFlux(mesh->Nlocal, (dfloat)0);
    platform->app->bc->o_usrwrk.copyFrom(zeroFlux);
  }

  const int Nfp = mesh->Nfp;
  const int Nq = mesh->Nq;

  int nSamplesRequested = 4096;
  platform->options.getArgs("RADIATION NSAMPLES", nSamplesRequested);
  int seed = 0;
  platform->options.getArgs("RADIATION SEED", seed);
  const bool writeMatrix = platform->options.compareArgs("RADIATION WRITE MATRIX", "TRUE");
  const bool useCache = platform->options.compareArgs("RADIATION CACHE", "TRUE");

  const std::string casename = platform->options.getArgs("CASENAME");
  std::string outputFile = casename + "_radiation";
  platform->options.getArgs("RADIATION OUTPUT FILE", outputFile);

  const std::string groupFile = outputFile + "_viewfactors_groups.csv";
  const std::string patchFile = outputFile + "_viewfactors_patches.bin";
  const std::string hashFile = outputFile + ".hash";

  const int nStrata = static_cast<int>(std::ceil(std::sqrt((double)nSamplesRequested)));
  const int nSamples = nStrata * nStrata;

  std::string fingerprint;
  int cacheHit = 0;
  if (rank == 0 && useCache) {
    std::ostringstream fp;
    fp << "P=" << P << " nSamples=" << nSamples << " seed=" << seed << " hasObstruction=" << hasObstruction
       << " radiating=" << radiatingStr << " obstruction=" << obstructionStr;
    for (int k = 0; k < global.K; ++k) {
      if (!global.isRadiating[k] && !global.isObstruction[k]) {
        continue;
      }
      const dfloat *c = &global.coords[(size_t)k * 3 * Nfp];
      for (int n = 0; n < 3 * Nfp; ++n) {
        fp << "," << c[n];
      }
    }
    fingerprint = SHA1::from_string(fp.str());

    std::ifstream hf(hashFile);
    if (hf.good()) {
      std::string cached;
      std::getline(hf, cached);
      if (cached == fingerprint) {
        cacheHit = 1;
      }
    }
  }
  MPI_Bcast(&cacheHit, 1, MPI_INT, 0, comm);

  if (cacheHit) {
    // The runtime radiosity coupling (Radiation::step) needs F in memory,
    // which the cache-hit path otherwise skips computing entirely. Reload it
    // from the cached patch-level binary matrix instead; if that file isn't
    // available (writeMatrix was off on the run that created the cache),
    // fall through and recompute rather than silently leaving F_state empty.
    int reloadOk = 0;
    if (rank == 0) {
      std::ifstream pf(patchFile, std::ios::binary);
      if (pf.good()) {
        int header[2];
        pf.read(reinterpret_cast<char *>(header), sizeof(header));
        if (pf.good() && header[0] == P && header[1] == P) {
          std::vector<double> Fcached((size_t)P * P);
          pf.read(reinterpret_cast<char *>(Fcached.data()), Fcached.size() * sizeof(double));
          if (pf.good()) {
            F_state.assign((size_t)P * P, (dfloat)0);
            for (size_t idx = 0; idx < Fcached.size(); ++idx) {
              F_state[idx] = (dfloat)Fcached[idx];
            }
            J_state.assign(P, (dfloat)0);
            reloadOk = 1;
          }
        }
      }
    }
    MPI_Bcast(&reloadOk, 1, MPI_INT, 0, comm);

    if (reloadOk) {
      if (rank == 0) {
        printf("Radiation: cache hit, %s already up to date, skipping Monte Carlo pass\n", groupFile.c_str());
      }
      setupCalled = true;
      return;
    }

    if (rank == 0) {
      printf("Radiation: cache hit but %s is unavailable (writeMatrix was likely off); "
            "recomputing view factors so the runtime radiosity coupling has a matrix to use\n",
            patchFile.c_str());
    }
  }

  // ---- radiating patch device arrays ----
  std::vector<dfloat> patchX((size_t)P * Nfp), patchY((size_t)P * Nfp), patchZ((size_t)P * Nfp);
  std::vector<dfloat> patchRefNx(P), patchRefNy(P), patchRefNz(P);
  std::vector<int> patchFaceKey(P);
  for (int p = 0; p < P; ++p) {
    const int k = radiatingIndices[p];
    for (int n = 0; n < Nfp; ++n) {
      patchX[(size_t)p * Nfp + n] = global.coords[(size_t)k * 3 * Nfp + 3 * n + 0];
      patchY[(size_t)p * Nfp + n] = global.coords[(size_t)k * 3 * Nfp + 3 * n + 1];
      patchZ[(size_t)p * Nfp + n] = global.coords[(size_t)k * 3 * Nfp + 3 * n + 2];
    }
    patchRefNx[p] = global.refNormal[3 * k + 0];
    patchRefNy[p] = global.refNormal[3 * k + 1];
    patchRefNz[p] = global.refNormal[3 * k + 2];
    patchFaceKey[p] = k;
  }

  deviceMemory<dfloat> o_patchX(patchX), o_patchY(patchY), o_patchZ(patchZ);
  deviceMemory<dfloat> o_patchRefNx(patchRefNx), o_patchRefNy(patchRefNy), o_patchRefNz(patchRefNz);
  deviceMemory<int> o_patchFaceKey(patchFaceKey);

  // ---- barycentric Lagrange weights for the mesh's own GLL nodes ----
  std::vector<dfloat> gllzHost(Nq), baryWHost(Nq);
  {
    std::vector<double> gllzD(Nq);
    for (int i = 0; i < Nq; ++i) {
      gllzD[i] = (double)mesh->gllz[i];
    }
    for (int j = 0; j < Nq; ++j) {
      double w = 1.0;
      for (int k = 0; k < Nq; ++k) {
        if (k != j) {
          w *= (gllzD[j] - gllzD[k]);
        }
      }
      baryWHost[j] = (dfloat)(1.0 / w);
      gllzHost[j] = (dfloat)gllzD[j];
    }
  }
  deviceMemory<dfloat> o_gllz(gllzHost), o_baryW(baryWHost);

  // ---- BVH over obstruction patches (host build, device traversal only) ----
  int nBvhNodes = 0;
  int bvhRoot = -1;
  deviceMemory<dfloat> o_bvhMin(1), o_bvhMax(1);
  deviceMemory<int> o_bvhLeft(1), o_bvhRight(1), o_bvhTriStart(1), o_bvhTriCount(1);
  deviceMemory<dfloat> o_triV0(1), o_triV1(1), o_triV2(1);
  deviceMemory<int> o_triPatchKey(1);

  if (hasObstruction && !obstructionIndices.empty()) {
    std::vector<std::vector<dfloat>> obsCoords(obstructionIndices.size());
    std::vector<int> obsKeys(obstructionIndices.size());
    for (size_t m = 0; m < obstructionIndices.size(); ++m) {
      const int k = obstructionIndices[m];
      obsCoords[m].assign(global.coords.begin() + (size_t)k * 3 * Nfp,
                          global.coords.begin() + (size_t)(k + 1) * 3 * Nfp);
      obsKeys[m] = k;
    }

    auto flat = RadiationBVH::build(obsCoords, obsKeys, Nq);
    nBvhNodes = static_cast<int>(flat.nodes.size());

    if (nBvhNodes > 0) {
      bvhRoot = nBvhNodes - 1;

      std::vector<dfloat> bvhMinH((size_t)nBvhNodes * 3), bvhMaxH((size_t)nBvhNodes * 3);
      std::vector<int> bvhLeftH(nBvhNodes), bvhRightH(nBvhNodes), bvhTriStartH(nBvhNodes), bvhTriCountH(nBvhNodes);
      for (int n = 0; n < nBvhNodes; ++n) {
        for (int d = 0; d < 3; ++d) {
          bvhMinH[3 * n + d] = flat.nodes[n].bmin[d];
          bvhMaxH[3 * n + d] = flat.nodes[n].bmax[d];
        }
        bvhLeftH[n] = flat.nodes[n].left;
        bvhRightH[n] = flat.nodes[n].right;
        bvhTriStartH[n] = flat.nodes[n].triStart;
        bvhTriCountH[n] = flat.nodes[n].triCount;
      }

      const int nTris = static_cast<int>(flat.triangles.size());
      std::vector<dfloat> triV0H((size_t)nTris * 3), triV1H((size_t)nTris * 3), triV2H((size_t)nTris * 3);
      std::vector<int> triPatchKeyH(nTris);
      for (int t = 0; t < nTris; ++t) {
        for (int d = 0; d < 3; ++d) {
          triV0H[3 * t + d] = flat.triangles[t].v0[d];
          triV1H[3 * t + d] = flat.triangles[t].v1[d];
          triV2H[3 * t + d] = flat.triangles[t].v2[d];
        }
        triPatchKeyH[t] = flat.triangles[t].patchKey;
      }

      o_bvhMin = deviceMemory<dfloat>(bvhMinH);
      o_bvhMax = deviceMemory<dfloat>(bvhMaxH);
      o_bvhLeft = deviceMemory<int>(bvhLeftH);
      o_bvhRight = deviceMemory<int>(bvhRightH);
      o_bvhTriStart = deviceMemory<int>(bvhTriStartH);
      o_bvhTriCount = deviceMemory<int>(bvhTriCountH);
      o_triV0 = deviceMemory<dfloat>(triV0H);
      o_triV1 = deviceMemory<dfloat>(triV1H);
      o_triV2 = deviceMemory<dfloat>(triV2H);
      o_triPatchKey = deviceMemory<int>(triPatchKeyH);
    }
  }

  // ---- partition the P*(P+1)/2 upper-triangle pair workload across ranks ----
  const long long totalPairs = (long long)P * (P + 1) / 2;

  int nRanks;
  MPI_Comm_size(comm, &nRanks);
  const long long pairsPerRank = totalPairs / nRanks;
  const long long remainder = totalPairs % nRanks;
  const long long myStart = (long long)rank * pairsPerRank + std::min<long long>(rank, remainder);
  const long long myCount = pairsPerRank + (rank < remainder ? 1 : 0);
  const long long myEnd = myStart + myCount;

  std::vector<int> localPairI, localPairJ;
  localPairI.reserve((size_t)myCount);
  localPairJ.reserve((size_t)myCount);
  {
    long long cursor = 0;
    for (int i = 0; i < P && cursor < myEnd; ++i) {
      const long long rowCount = P - i;
      const long long rowStart = cursor;
      const long long rowEnd = cursor + rowCount;
      if (rowEnd > myStart && rowStart < myEnd) {
        const long long overlapStart = std::max(rowStart, myStart);
        const long long overlapEnd = std::min(rowEnd, myEnd);
        for (long long p = overlapStart; p < overlapEnd; ++p) {
          const int j = i + static_cast<int>(p - rowStart);
          localPairI.push_back(i);
          localPairJ.push_back(j);
        }
      }
      cursor = rowEnd;
    }
  }

  const dlong nPairsLocal = static_cast<dlong>(localPairI.size());

  deviceMemory<int> o_pairI(nPairsLocal ? localPairI : std::vector<int>(1, 0));
  deviceMemory<int> o_pairJ(nPairsLocal ? localPairJ : std::vector<int>(1, 0));

  deviceMemory<double> o_rawIJLocal(std::max<dlong>(nPairsLocal, 1));
  deviceMemory<double> o_areaEstLocal(P);
  {
    std::vector<double> zerosP(P, 0.0);
    o_areaEstLocal.copyFrom(zerosP);
  }

  if (rank == 0) {
    printf("Radiation: sampling %lld view-factor pairs with %d MC samples/pair...\n", totalPairs, nSamples);
  }

  if (nPairsLocal > 0) {
    samplePairsKernel(nPairsLocal,
                      (dlong)nSamples,
                      (dlong)nStrata,
                      (unsigned int)seed,
                      (dlong)nBvhNodes,
                      bvhRoot,
                      o_pairI,
                      o_pairJ,
                      o_patchX,
                      o_patchY,
                      o_patchZ,
                      o_patchRefNx,
                      o_patchRefNy,
                      o_patchRefNz,
                      o_patchFaceKey,
                      o_gllz,
                      o_baryW,
                      o_bvhMin,
                      o_bvhMax,
                      o_bvhLeft,
                      o_bvhRight,
                      o_bvhTriStart,
                      o_bvhTriCount,
                      o_triV0,
                      o_triV1,
                      o_triV2,
                      o_triPatchKey,
                      o_rawIJLocal,
                      o_areaEstLocal);
  }

  // ---- gather results and combine across ranks ----
  std::vector<double> rawIJGlobal((size_t)totalPairs, 0.0);
  if (nPairsLocal > 0) {
    std::vector<double> rawIJLocalHost(nPairsLocal);
    o_rawIJLocal.copyTo(rawIJLocalHost);
    std::copy(rawIJLocalHost.begin(), rawIJLocalHost.end(), rawIJGlobal.begin() + myStart);
  }

  std::vector<double> areaEstGlobal(P, 0.0);
  o_areaEstLocal.copyTo(areaEstGlobal);

  MPI_Allreduce(MPI_IN_PLACE, rawIJGlobal.data(), (int)totalPairs, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, areaEstGlobal.data(), P, MPI_DOUBLE, MPI_SUM, comm);

  auto packedIndex = [P](int i, int j) -> long long {
    if (i > j) {
      std::swap(i, j);
    }
    return (long long)i * P - (long long)i * (i - 1) / 2 + (j - i);
  };

  if (rank == 0) {
    std::vector<double> F((size_t)P * P, 0.0);
    for (int i = 0; i < P; ++i) {
      for (int j = i; j < P; ++j) {
        const double raw = rawIJGlobal[packedIndex(i, j)];
        const double Ai = areaEstGlobal[i];
        const double Aj = areaEstGlobal[j];
        F[(size_t)i * P + j] = (Ai > 0) ? raw / Ai : 0.0;
        F[(size_t)j * P + i] = (Aj > 0) ? raw / Aj : 0.0;
      }
    }

    // Retained in memory (independent of writeMatrix, which only controls
    // the on-disk copy) for Radiation::step()'s runtime radiosity coupling.
    F_state.assign((size_t)P * P, (dfloat)0);
    for (size_t idx = 0; idx < F.size(); ++idx) {
      F_state[idx] = (dfloat)F[idx];
    }
    J_state.assign(P, (dfloat)0);

    std::vector<double> groupArea(nGroups, 0.0);
    for (int p = 0; p < P; ++p) {
      const int gi = patchGroupIdx[p];
      groupArea[gi] += areaEstGlobal[p];
    }

    std::vector<double> groupFlux((size_t)nGroups * nGroups, 0.0);
    for (int i = 0; i < P; ++i) {
      const int gi = patchGroupIdx[i];
      for (int j = 0; j < P; ++j) {
        const int gj = patchGroupIdx[j];
        groupFlux[(size_t)gi * nGroups + gj] += areaEstGlobal[i] * F[(size_t)i * P + j];
      }
    }

    std::ofstream gf(groupFile);
    gf << std::setprecision(15);
    gf << "g,h,F_gh,A_g\n";
    for (int gi = 0; gi < nGroups; ++gi) {
      for (int gj = 0; gj < nGroups; ++gj) {
        const double Fgh = (groupArea[gi] > 0) ? groupFlux[(size_t)gi * nGroups + gj] / groupArea[gi] : 0.0;
        gf << groups[gi] << "," << groups[gj] << "," << Fgh << "," << groupArea[gi] << "\n";
      }
    }
    gf.close();

    if (writeMatrix) {
      std::ofstream pf(patchFile, std::ios::binary);
      int header[2] = {P, P};
      pf.write(reinterpret_cast<const char *>(header), sizeof(header));
      pf.write(reinterpret_cast<const char *>(F.data()), F.size() * sizeof(double));
      pf.close();
    }

    if (useCache) {
      std::ofstream hf(hashFile, std::ios::trunc);
      hf << fingerprint;
      hf.close();
    }

    std::string msg = groupFile;
    if (writeMatrix) {
      msg += " and " + patchFile;
    }
    printf("Radiation: wrote %s\n", msg.c_str());
  }

  setupCalled = true;
}

// Couples the (fixed) view-factor matrix F_state to the (time-varying)
// temperature field: gathers per-patch average temperature, solves the
// gray-diffuse radiosity system (warm-started from the previous call) for
// the net radiative flux leaving each patch, and scatters the result into
// platform->app->bc->o_usrwrk for the case's udfNeumann to consume. No-op
// except every RADIATION updateFrequency steps -- the flux from the last
// update is left in place in between.
void Radiation::step(double time, int tstep)
{
  if (!setupCalled) {
    return;
  }
  if (updateFrequency_state <= 0 || tstep % updateFrequency_state != 0) {
    return;
  }

  MPI_Comm comm = platform->comm.mpiComm();
  const int rank = platform->comm.mpiRank();
  auto mesh = nrs->meshV;
  const int P = P_state;

  nekrsCheck(!nrs->scalar || nrs->scalar->nameToIndex.find("temperature") == nrs->scalar->nameToIndex.end(),
            comm,
            EXIT_FAILURE,
            "%s\n",
            "[RADIATION] step() requires a [TEMPERATURE] scalar field to couple to.");

  // ---- this rank's per-local-patch average temperature ----
  auto o_T = nrs->scalar->o_solution("temperature");
  std::vector<dfloat> Thost(mesh->Nlocal);
  o_T.copyTo(Thost, mesh->Nlocal);

  const int nLocalRad = static_cast<int>(localRadiatingPatches_state.size());
  std::vector<int> localP(std::max(nLocalRad, 1));
  std::vector<dfloat> localT(std::max(nLocalRad, 1));
  for (int i = 0; i < nLocalRad; ++i) {
    const auto &patch = localRadiatingPatches_state[i];
    dfloat sum = 0;
    for (dlong idv : patch.idxVol) {
      sum += Thost[idv];
    }
    localP[i] = patch.p;
    localT[i] = sum / (dfloat)patch.idxVol.size();
  }

  // ---- rank 0 collects (patch index, T) pairs from every rank ----
  int nRanks;
  MPI_Comm_size(comm, &nRanks);
  std::vector<int> counts(nRanks), displs(nRanks, 0);
  MPI_Gather(&nLocalRad, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  if (rank == 0) {
    int total = 0;
    for (int r = 0; r < nRanks; ++r) {
      displs[r] = total;
      total += counts[r];
    }
  }

  std::vector<int> allP(rank == 0 ? P : 0);
  std::vector<dfloat> allT(rank == 0 ? P : 0);
  MPI_Gatherv(localP.data(), nLocalRad, MPI_INT, allP.data(), counts.data(), displs.data(), MPI_INT, 0, comm);
  MPI_Gatherv(localT.data(), nLocalRad, MPI_DFLOAT, allT.data(), counts.data(), displs.data(), MPI_DFLOAT, 0, comm);

  // ---- rank 0: gray-diffuse radiosity solve, warm-started from J_state ----
  std::vector<dfloat> q(P, (dfloat)0);
  if (rank == 0) {
    std::vector<dfloat> T(P, (dfloat)0);
    for (int i = 0; i < P; ++i) {
      T[allP[i]] = allT[i];
    }

    std::vector<dfloat> Eb(P);
    for (int i = 0; i < P; ++i) {
      const dfloat Ti2 = T[i] * T[i];
      Eb[i] = stefanBoltzmann_state * Ti2 * Ti2;
    }

    for (int iter = 0; iter < radiosityMaxIters_state; ++iter) {
      dfloat maxDelta = 0;
      for (int i = 0; i < P; ++i) {
        dfloat sum = 0;
        const dfloat *Fi = &F_state[(size_t)i * P];
        for (int j = 0; j < P; ++j) {
          sum += Fi[j] * J_state[j];
        }
        const dfloat eps = patchEmissivity_state[i];
        const dfloat Jnew = eps * Eb[i] + (1 - eps) * sum;
        maxDelta = std::max(maxDelta, (dfloat)std::fabs(Jnew - J_state[i]));
        J_state[i] = Jnew;
      }
      if (maxDelta < radiosityTolerance_state) {
        break;
      }
    }

    // Net radiative flux leaving patch i, W/m^2 (positive: patch i is a net
    // radiator, e.g. hotter than what it sees).
    for (int i = 0; i < P; ++i) {
      dfloat sum = 0;
      const dfloat *Fi = &F_state[(size_t)i * P];
      for (int j = 0; j < P; ++j) {
        sum += Fi[j] * J_state[j];
      }
      q[i] = patchEmissivity_state[i] * (Eb[i] - sum);
    }

    const int nGroups = static_cast<int>(groups_state.size());
    std::vector<dfloat> groupQSum(nGroups, 0);
    std::vector<dfloat> groupTSum(nGroups, 0);
    std::vector<int> groupCount(nGroups, 0);
    for (int i = 0; i < P; ++i) {
      const int gi = patchGroupIdx_state[i];
      groupQSum[gi] += q[i];
      groupTSum[gi] += T[i];
      groupCount[gi]++;
    }
    printf("Radiation: t=%g step=%d net radiative flux / avg temperature by boundary-ID group (W/m^2, K, avg over "
          "patches):\n",
          time,
          tstep);
    for (int gi = 0; gi < nGroups; ++gi) {
      printf("  bID=%d: %g  %g\n",
            groups_state[gi],
            groupCount[gi] ? groupQSum[gi] / groupCount[gi] : (dfloat)0,
            groupCount[gi] ? groupTSum[gi] / groupCount[gi] : (dfloat)0);
    }
  }

  MPI_Bcast(q.data(), P, MPI_DFLOAT, 0, comm);

  // ---- scatter q back into bc->o_usrwrk at each local patch's nodes ----
  std::vector<dfloat> fluxHost(mesh->Nlocal, (dfloat)0);
  for (const auto &patch : localRadiatingPatches_state) {
    for (dlong idv : patch.idxVol) {
      fluxHost[idv] = q[patch.p];
    }
  }
  platform->app->bc->o_usrwrk.copyFrom(fluxHost);
}
