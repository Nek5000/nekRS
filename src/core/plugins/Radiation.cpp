#include "nrs.hpp"
#include "platform.hpp"
#include "Radiation.hpp"
#include "RadiationBVH.hpp"
#include "mesh3D.h"
#include "sha1.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>

namespace
{

nrs_t *nrs;

bool buildKernelCalled = false;
bool setupCalled = false;

occa::kernel samplePairsKernel;

struct RadiationPatch {
  std::vector<dfloat> coords; // 3*Nfp, laid out coords[3*n+0/1/2]
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

      for (int n = 0; n < mesh->Nfp; ++n) {
        const dlong idM = mesh->vmapM[e * mesh->Nfaces * mesh->Nfp + f * mesh->Nfp + n];
        patch.coords[3 * n + 0] = x[idM];
        patch.coords[3 * n + 1] = y[idM];
        patch.coords[3 * n + 2] = z[idM];
      }

      const dlong sid = e * mesh->Nfaces * mesh->Nfp + f * mesh->Nfp + 0;
      mesh->o_sgeo.copyTo(sgeoNode.data(), mesh->Nsgeo, sid * mesh->Nsgeo);
      patch.refNx = sgeoNode[NXID];
      patch.refNy = sgeoNode[NYID];
      patch.refNz = sgeoNode[NZID];

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
                               MPI_Comm comm)
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

  int nRanks;
  MPI_Comm_size(comm, &nRanks);

  std::vector<int> counts(nRanks);
  MPI_Allgather(&nLocal, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);

  std::vector<int> displs(nRanks);
  int total = 0;
  for (int r = 0; r < nRanks; ++r) {
    displs[r] = total;
    total += counts[r];
  }

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
  auto mesh = nrs->fluid->mesh;

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

  const auto global = allgatherPatches(localPatches, isRadiatingFlag, isObstructionFlag, mesh->Nfp, comm);

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
    if (rank == 0) {
      printf("Radiation: cache hit, %s already up to date, skipping Monte Carlo pass\n", groupFile.c_str());
    }
    setupCalled = true;
    return;
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

    std::set<int> groupIDSet;
    for (int p = 0; p < P; ++p) {
      groupIDSet.insert(global.boundaryID[radiatingIndices[p]]);
    }
    const std::vector<int> groups(groupIDSet.begin(), groupIDSet.end());
    const int nGroups = static_cast<int>(groups.size());

    std::vector<int> patchGroupIdx(P);
    std::vector<double> groupArea(nGroups, 0.0);
    for (int p = 0; p < P; ++p) {
      const int bID = global.boundaryID[radiatingIndices[p]];
      const int gi = static_cast<int>(std::lower_bound(groups.begin(), groups.end(), bID) - groups.begin());
      patchGroupIdx[p] = gi;
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
