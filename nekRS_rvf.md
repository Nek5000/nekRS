# nekRS Radiative View-Factor Module — Session Log

This file is a written record of the chat session that designed and implemented
a Monte Carlo radiative view-factor module for nekRS, ported from the algorithm
in `RadiativeViewFactor.jl`.

---

## User

I would like to build a radiative heat transfer module for nekRS: i.e. something
that you can turn on in the .par file for a case that computes view factors for
specified mesh elements/groups like RadiativeViewFactor.jl.
It should exclusively use monte-carlo ray tracing methods to calculate the view
factor, and all the code should be written in OCCA/C++ so that it integrates
naturally with nekRS.
After the view factor code is written and integrated, we'll talk about how to
implement radiative heat transfer during a CFD run.

## Assistant

Launched two research passes in parallel before designing anything:

1. **RadiativeViewFactor.jl algorithm** — confirmed it is **area-to-area
   stratified Monte Carlo sampling**, not hemisphere ray casting from a point:
   - For patch pair (i,j): draw N point pairs, one on each patch's reference
     domain, on a stratified (jittered) `⌊√N⌋×⌊√N⌋` grid.
   - Map each reference point to physical space via the patch's own
     isoparametric shape functions (curved, not flat) to get position, normal,
     and differential area.
   - Kernel `K = cosθᵢ·cosθⱼ / (π r²)`, zero if either cosine ≤ 0.
   - Occlusion: segment visibility test against a BVH built from a
     user-specified obstruction set; source/destination patches self-excluded;
     skipped entirely if no obstruction set configured.
   - `F_ij = raw_ij / A_i`, with `A_i` from patch i's own dedicated self-pair
     (`i==i`) MC run — this makes `raw_ij` symmetric by construction, so
     reciprocity is exact, not enforced post hoc.
   - GPU version: one thread(-group) per pair, GPU-safe inline PRNG, 32-bit
     preferred (Metal/Apple GPU compatibility).

2. **nekRS plugin architecture** — confirmed:
   - Optional physics modules live in `src/core/plugins/` (generic, e.g.
     `tavg`, `lpm`) or `src/app/nrs/plugins/` (NS-specific, e.g. `RANSktau`,
     `lowMach`), each as a `.hpp`/`.cpp` pair plus a sibling `kernels/*.okl`.
   - `.par` parsing goes through `inipp`, with per-section parsers dispatched
     from `src/platform/par/par.cpp`, which enforces a **hard whitelist** of
     valid sections/keys (`validSections`, `getValidKeys()`).
   - OCCA kernels are registered/built via `platform->kernelRequests`
     (two-phase: `add` during a "REGISTER ONLY" pass, then `load`). `.okl`
     files need no CMake registration — auto-installed by
     `cmake/install_kernels.cmake` globbing.
   - `UDF_AutoLoadPlugins` is auto-generated boilerplate (`udfMake.hpp`)
     keyed on each plugin header's include guard — adding an entry there
     means a case's `.udf` only needs `#include "Radiation.hpp"` to get
     kernel registration for free. The actual `setup()` call still needs one
     line in `UDF_Setup()`, matching every other optional module's
     convention.
   - Boundary faces: `mesh_t::EToB` (element-to-boundary-ID), `vmapM`
     (face-node → volume-node index), `Nfp = Nq*Nq` nodes/face (hex-only
     mesh). `o_sgeo`'s `NXID/NYID/NZID/WSJID` are defined **only at the fixed
     GLL nodes**, not usable for arbitrary MC sample points.
   - No existing radiation/view-factor code anywhere in the repo.
   - nekRS is MPI-distributed; boundary faces are partitioned across ranks.

Then entered plan mode and asked one clarifying question before designing:

**Q: face geometry fidelity for the MC sampling — curvature-accurate or a
flat/bilinear shortcut for v1?**
→ User chose **curvature-accurate** (map samples through the true polynomial
face geometry via tensor-product Lagrange interpolation over the face's own
GLL nodes, generalizing the Julia package's Quad8 isoparametric approach to
arbitrary order).

Spawned a Plan agent to turn this into a concrete implementation plan, then
verified its two most load-bearing claims directly against source:
- `par.cpp:434` really does hard-reject unregistered sections — confirmed.
- `FaceNodesHex3D` (`meshBasisHex3D.cpp:55-93`) builds each face's node list
  by scanning the volume's own tensor-product node order (`n = i + j*Nq +
  k*Nq*Nq`, confirmed in `NodesHex3D`), so every face's `Nfp` list is a valid
  row-major 2D grid over its two free reference directions — exactly what the
  Lagrange-mapping kernel needs. Confirmed for all 6 faces, not just
  axis-aligned ones.

Wrote the plan to `/Users/coxea3/.claude/plans/jaunty-frolicking-kazoo.md` and
called `ExitPlanMode`.

## User (plan-mode rejection)

> The plan looks good, but instead of referring to sections/files as
> 'viewfactor' or similar, name things 'radiation' or similar

## Assistant

Renamed every identifier throughout the plan: `[VIEWFACTOR]` → `[RADIATION]`,
`ViewFactor.hpp/cpp` → `Radiation.hpp/cpp`, `ViewFactorBVH` → `RadiationBVH`,
`ViewFactor::setup` → `Radiation::setup`, output filenames, example directory
name, etc. — with the rationale that view factors are one component of a
broader radiation module this will grow into (heat-transfer coupling is the
explicitly-planned next phase), so the name shouldn't be view-factor-specific.
Re-submitted via `ExitPlanMode`; the user approved.

### Approved plan (final)

**Context.** RadiativeViewFactor.jl already computes view factors via
stratified area-to-area Monte Carlo sampling with BVH occlusion. Phase 1 here
is porting that method into nekRS itself as a `.par`-toggled module operating
directly on nekRS's own (possibly curved, high-order) mesh — no separate
meshing/export step. Only view-factor computation and output are in scope;
coupling into an actual radiative heat-transfer boundary condition is a later,
separate effort.

**Key structural facts** (verified, not assumed): hex-only `mesh_t`; a
boundary face is `(e,f)` with `EToB[f+e*Nfaces]>0`; its `Nfp` nodal coords are
`vmapM[e*Nfaces*Nfp+f*Nfp+n]` into `o_x/o_y/o_z`; each face's node list is a
valid row-major tensor grid over its two free reference directions (verified
against `meshBasisHex3D.cpp`); `.par` has a hard section/key whitelist;
`UDF_AutoLoadPlugins` boilerplate is auto-generated from a header-include-guard
table; nekRS is MPI-distributed so patch geometry must be globally replicated
(`MPI_Allgatherv`) and the O(P²) pair workload striped across ranks
(`MPI_Allreduce` to combine); no existing helper computes barycentric Lagrange
weights for arbitrary-point evaluation (new, ~20 lines, standard formula).

**Design:**

1. **`.par` section `[RADIATION]`** — `radiatingBoundaryIDs` (required, also
   the aggregation group unit), `obstructionBoundaryIDs` (optional; no
   occlusion test at all if absent), `nSamples`, `seed`, `writeMatrix`,
   `outputFile`, `cache`. Registered in `par.cpp`'s validator tables via a new
   `parseRadiation.hpp`, modeled on `parseNeknek.hpp`.

2. **Patch representation** — one `(e,f)` boundary face = one patch, for both
   radiating and obstruction sets (a face may be in both). Local faces
   gathered per-rank, then globally replicated via `MPI_Allgatherv`; warn/hard
   error past a configurable patch-count threshold since this targets a
   "surface enclosure" subset, not the whole boundary.

3. **BVH (occlusion structure)** — triangulate each obstruction patch's
   `Nq×Nq` nodal grid into `2(Nq-1)²` flat triangles (curvature resolution
   scales with mesh order); build a flat/linear BVH **on the host**, once, at
   setup; device does traversal only. Skipped entirely if no obstruction set
   configured.

4. **OKL kernel** — one `@outer` group per pair `(i,j)` (upper triangle only,
   exploiting symmetry), matching the `@shared`/`@inner(0,1)` idiom already
   used in `RANSktau.okl` (cooperative load of both patches' nodal grids into
   shared memory, then `@inner` threads stream a partition of the N samples
   with per-thread accumulation + in-block reduction — deliberately not a
   literal one-thread-per-pair design, to avoid register pressure at
   realistic `Nq`/`N`). Per-sample device routine does a two-pass separable
   barycentric Lagrange contraction (Berrut/Trefethen form) over the shared
   nodal grid — O(Nq²) unavoidable minimum — producing position and, via the
   two tangent partials, an un-normalized normal whose magnitude is `dA_ref`
   directly (no extra Jacobian bookkeeping needed, since both the GLL node
   range and the MC reference domain are `[-1,1]²` by construction). RNG: a
   32-bit Wang-hash-style seed mix (not a literal 64-bit port — no `.okl` file
   in this codebase uses 64-bit int ops, a real portability risk on some
   OpenCL backends). Occlusion: stackless/small-fixed-stack BVH traversal,
   Möller–Trumbore clipped to the open segment so patch endpoints never
   self-occlude. Self-pairs (`i==j`) draw an independent second sample pair
   and are the sole writer of the area estimate — cross pairs never
   re-estimate area, so `F_ij = raw_ij/areaEst[i]` is well-defined and
   reciprocity-consistent.

5. **Host post-processing** — combine MPI-partitioned results
   (`MPI_Allreduce`), compute `F[i][j]` and its mirror, aggregate to
   boundary-ID group level (`A_g·F_{g→h} = Σ_{i∈g}Σ_{j∈h} A_i·F_ij`). Rank 0
   writes a group-level CSV (always) and a dense patch matrix (if
   `writeMatrix`).

6. **Caching** — SHA1 fingerprint over mesh geometry + config, stored as a
   sibling hash file; skips the whole MC pass on an unchanged case.

7. **File list** — `src/core/plugins/Radiation.{hpp,cpp}`,
   `RadiationBVH.{hpp,cpp}`, `kernels/Radiation.okl`,
   `src/platform/par/parseRadiation.hpp` (+ `par.cpp` registration),
   `udfMake.hpp` pluginTable entry, `cmake/core.cmake` source-list addition,
   plus a planned (not yet written) `examples/radiationPlates` validation
   case.

**Verification plan** — closed-form parallel/perpendicular plate view factors
(Modest, *Radiative Heat Transfer*) with `1/√N` MC convergence and
polynomial-order convergence checks on a curved plate; reciprocity check
(`A_i F_ij == A_j F_ji`); closure check on a closed enclosure
(`Σ_h F_gh == 1`); occlusion check (obstructed plates → `F_12` drops, → 0 when
fully blocked); determinism/cache-hit check.

## Assistant — implementation

Read the concrete reference files needed to match nekRS's real APIs exactly
(`RANSktau.hpp/.cpp/.okl`, `parseNeknek.hpp`, `par.cpp`, `mesh.h`,
`meshConnectFaceNodes3D.cpp`, `deviceMemory.hpp`, `udfMake.hpp`,
`cmake/core.cmake`, `platform.hpp`, `setupAide.hpp`, `sha1.hpp`,
`meshSurface.cpp`, `nrs.hpp`, `fluidSolver.hpp`, `kernelManager.hpp`,
`nekrsSys.hpp.in`, `meshBasisHex3D.cpp`), and confirmed via
`cmake/install_kernels.cmake` the exact installed-kernel-path convention
(`.../kernels/` path component stripped) so `buildKernel()`'s file path would
be correct.

Noted the nekRS repo was already on a branch named `rht` (radiative heat
transfer) with only trivial pre-existing noise (an editor `.gitignore` entry,
an unrelated submodule whitespace diff) — no prior radiation code, so this was
a clean build on the intended branch.

Wrote the module (see file list below), then wired it into the build.

### Files created

| File | Purpose |
|---|---|
| `src/platform/par/parseRadiation.hpp` | `[RADIATION]` `.par` section parser |
| `src/core/plugins/Radiation.hpp` | public API: `buildKernel`, `setup` |
| `src/core/plugins/Radiation.cpp` | host orchestration: parses config, gathers boundary faces into patches, MPI-replicates geometry (`MPI_Allgatherv`), computes barycentric Lagrange weights, partitions the O(P²) pair workload across ranks, launches the kernel, combines results (`MPI_Allreduce`), aggregates to boundary-ID groups, writes CSV/binary output, SHA1-based caching |
| `src/core/plugins/RadiationBVH.hpp` | BVH data structures + `build()` declaration |
| `src/core/plugins/RadiationBVH.cpp` | host-side triangulation of curved obstruction patches + median-split BVH build |
| `src/core/plugins/kernels/Radiation.okl` | the MC kernel: stratified sampling, curvature-accurate tensor-product barycentric Lagrange face mapping, stackless BVH/Möller–Trumbore occlusion test, counter-based 32-bit RNG |

### Files edited

| File | Change |
|---|---|
| `src/platform/par/par.cpp` | added `radiationKeys`, registered `"radiation"` in `validSections`, wired into `makeStringsLowerCase()`/`getValidKeys()`, `#include "parseRadiation.hpp"`, call site in `Par::parse` |
| `src/core/udf/udfMake.hpp` | added `{"nekrs_Radiation_hpp_", "Radiation::buildKernel"}` to the autoload `pluginTable` |
| `cmake/core.cmake` | added `Radiation.cpp` and `RadiationBVH.cpp` to `CORE_SOURCES` |

### Deliberate deviations from the plan, made during implementation

- **Lowered the patch-count cap from ~50k to 8000.** The plan's figure didn't
  account for O(P²) dense-matrix memory — 50k patches would need ~10GB just
  for the raw pair matrix; 8000 keeps it under ~256MB.
- **Explicit MPI pair-workload striping.** Each rank owns a contiguous slice
  of the upper-triangle pair enumeration (computed via row-cumulative-count
  inversion, avoiding ever materializing the full O(P²) pair list), combined
  via `MPI_Allreduce`. This wasn't fully resolved in the plan's kernel-launch
  section and was worked out during implementation.
- **Reference-outward-normal sign correction.** Not explicitly called out in
  the plan: since the Lagrange-mapping's cross-product normal direction
  depends on an arbitrary internal parametrization, each patch's true outward
  direction (read once from `mesh->o_sgeo`'s `NXID/NYID/NZID` at one face
  node — the only reliable outward-consistent source nekRS provides) is used
  to sign-correct the kernel's computed normal at every sample point. Without
  this, the cosine-sign convention in the MC kernel would be wrong for
  faces whose internal tensor-grid orientation happens to give an
  inward-pointing cross product.
- **Stateless counter-based RNG instead of a literal xorshift32 stream.** Each
  draw is `hash(seed, patchI, patchJ, drawSlot)` rather than a mutable
  per-thread xorshift32 state advanced call-by-call — simpler to get right
  inside nested OKL conditionals, equally GPU-safe/stateless, same intent as
  the plan's RNG design.

### Not done / not verified

- **No compile attempted.** `nekrsSys.hpp` is templated (`.in`) and only
  generated by CMake configure; a full build (OCCA submodule + MPI + nekRS)
  is a substantial, environment-sensitive process not attempted in this
  session without explicit direction to do so. Everything was checked
  statically against real usages elsewhere in the codebase (types, argument
  order across the `.cpp`↔`.okl` boundary, API signatures), but a real build
  pass is the next step.
- **No validation example case** (`examples/radiationPlates`). Generating a
  mesh needs `genbox`/Nek5000 tooling not confirmed working in this
  environment, so no untested case was fabricated.
- **Portability caveat** (flagged in the `.okl` comments): the helper device
  functions take `@shared`-memory pointers as plain arguments, which works on
  CUDA/HIP/Serial but is a known rough edge on strict OpenCL backends.

### Open question for the user

Attempt a build now and iterate on compiler errors, or set up the validation
case (parallel-plates closed-form check) first?

## User

Please write out this whole chat to a file nekRS_rvf.md in the nekRS folder

## Assistant

Wrote this file.
