# Radiative View-Factor Module — Code Changes

Full list of files added or modified in this repo (branch `rht`) to implement
the Monte Carlo radiative view-factor module (`[RADIATION]` in `.par`).
Validated against the closed-form Hottel/Feingold parallel-plate formula
(0.004% error) via `examples/radiationPlates`.

---

## New files

### `src/core/plugins/Radiation.hpp`
Public API: `Radiation::buildKernel(occa::properties)` (kernel registration,
called automatically via the UDF autoload table) and `Radiation::setup()`
(the one line a case's `UDF_Setup()` calls to run the module).

### `src/core/plugins/Radiation.cpp`
Host-side orchestration — the bulk of the module's logic:
- Parses `[RADIATION]` options (`radiatingBoundaryIDs`, `obstructionBoundaryIDs`,
  `nSamples`, `seed`, `writeMatrix`, `outputFile`, `cache`).
- Walks `mesh->EToB` to find local boundary faces matching the configured
  boundary IDs, gathers each face's `Nfp` nodal coordinates (via `vmapM`) and
  a reference outward-into-the-enclosure normal (read from `mesh->o_sgeo`,
  **negated** — sgeo's normal points out of the CFD domain, but a radiating
  surface's normal must point into the enclosure it bounds).
- `MPI_Allgatherv`s every rank's local patches into an identical global patch
  list on every rank (boundary faces are MPI-partitioned in nekRS), giving a
  deterministic global index used as a face-identity key.
- Computes barycentric Lagrange weights for the mesh's own GLL nodes (new
  code — no existing nekRS helper does this).
- Builds the BVH over obstruction-boundary patches (if configured) via
  `RadiationBVH::build`.
- Partitions the `P·(P+1)/2` upper-triangle patch-pair workload across MPI
  ranks (each rank computes a contiguous slice; no full O(P²) list is ever
  materialized), launches `radiationSamplePairs`, and `MPI_Allreduce`s the
  per-rank partial results into the full matrix on every rank.
- Post-processes into `F[i][j]` and an area-weighted boundary-ID-group
  aggregation, writes `<case>_radiation_viewfactors_groups.csv` (always) and
  `<case>_radiation_viewfactors_patches.bin` (if `writeMatrix`), with a
  SHA1-fingerprint cache (`<case>_radiation.hash`) that skips recomputation
  on an unchanged mesh/config.
- Uses `nrs->meshV`, not `nrs->fluid->mesh` — makes the module work whether
  or not `[FLUID VELOCITY]` is configured (a case with `FLUID=FALSE` leaves
  `nrs->fluid` null).

### `src/core/plugins/RadiationBVH.hpp` / `.cpp`
Host-side, one-time BVH build over the obstruction-patch triangle soup:
triangulates each obstruction patch's `Nq×Nq` curved nodal grid into
`2(Nq-1)²` flat triangles, then builds a median-split linear BVH (simple
top-down recursive build, leaf size 4) flattened into arrays for device
traversal. Kept in its own translation unit for isolated testability.

### `src/core/plugins/kernels/Radiation.okl`
The OCCA/OKL Monte Carlo kernel — `radiationSamplePairs`, one `@outer` group
per radiating-patch pair `(i,j)` (upper triangle incl. self-pairs):
- Cooperatively loads both patches' `Nq×Nq` nodal grids into `@shared`
  memory once per pair, then `@inner` threads each stream a strided subset
  of the stratified MC samples.
- `evalFaceLagrange`: maps a `[-1,1]²` reference point through the face's
  true (possibly curved) geometry via tensor-product barycentric Lagrange
  interpolation (Berrut/Trefethen form) over the shared nodal grid, giving
  position, normal, and differential area in one pass.
- `radiationOccluded`: stackless (fixed local-array stack) BVH traversal +
  Möller–Trumbore segment/triangle test for the optional occlusion check,
  self-excluding triangles belonging to either endpoint patch.
- `radiationRand`/`radiationWangHash`: a stateless, counter-based 32-bit
  hash RNG keyed on `(seed, patchI, patchJ, slot)` — no mutable per-thread
  RNG state, safe across the `@outer`/`@inner` structure.
- Self-pairs (`i==j`) draw an independent second sample set and are the sole
  writer of that patch's MC-estimated area; cross pairs only accumulate the
  raw kernel integral.

### `src/platform/par/parseRadiation.hpp`
`.par` section parser for `[RADIATION]`, modeled on `parseNeknek.hpp`.
Populates `platform->options` (`RADIATION ...` keys) from the section.

### `examples/radiationPlates/`
New validation case — two parallel unit-square plates, unit gap, no flow
solver:
- `radiationPlates.box` — Nek5000 `genbox` input (6×6×1 element box mesh,
  `[0,1]³`; boundary labels `O` for the 4 side walls (excluded), `W`/`v` for
  the z=0/z=1 plates, each getting its own boundary ID).
- `radiationPlates.re2` — the generated binary mesh.
- `radiationPlates.par` — `[GENERAL]` + `[MESH]` + `[RADIATION]` only (no
  `[FLUID VELOCITY]` needed, confirming the `meshV` fix above).
- `radiationPlates.udf` — the minimal 2-line integration:
  `#include "Radiation.hpp"` + `Radiation::setup()` in `UDF_Setup()`.
- `probe.par` / `probe.udf` — a throwaway diagnostic case (dumps `EToB` +
  face centroids) used to empirically determine which boundary ID ends up
  on which geometric face, since that mapping isn't independently
  derivable from `genbox`'s documentation. Kept for reference; not part of
  the module itself.

Result: `F_1→2 = 0.199833` vs. the closed-form value `0.199825` (0.004%
error); `F_1→2 = F_2→1` to 15 digits (reciprocity, by construction); `F_1,1`
and `F_2,2` both numerically zero (correct — flat plates can't see
themselves).

---

## Modified files

### `src/platform/par/par.cpp`
Registers `[RADIATION]` as a valid `.par` section: adds `radiationKeys` (the
7 accepted keys), adds `"radiation"` to `validSections`, wires it into
`makeStringsLowerCase()` and `getValidKeys()`, includes `parseRadiation.hpp`,
and calls `parseRadiationSection(...)` from `Par::parse`. Without this, `.par`'s
hard section/key whitelist rejects `[RADIATION]` outright.

### `src/core/udf/udfMake.hpp`
Adds `{"nekrs_Radiation_hpp_", "Radiation::buildKernel"}` to the
autogenerated `UDF_AutoLoadPlugins` table, so a case's `.udf` merely
`#include`ing `"Radiation.hpp"` is enough to get kernel registration for
free (matching the existing `RANSktau`/`lowMach`/`tavg` convention).

### `cmake/core.cmake`
Adds `src/core/plugins/Radiation.cpp` and `RadiationBVH.cpp` to
`CORE_SOURCES` so they're compiled into `libnekrs`/`libnekrs-fp32`.
`Radiation.okl` needed no CMake registration — `.okl` kernel files are
auto-discovered and installed by `cmake/install_kernels.cmake`.

---

## Bug fixes found via the validation case

These were caught only by actually running the module end-to-end and
checking the result against a closed-form value — worth flagging since
they're the kind of thing that would silently produce wrong physics
without this validation:

1. **Inverted normal convention** (`Radiation.cpp`): `mesh->o_sgeo`'s normal
   points out of the CFD domain; a radiating surface needs its normal
   pointing into the enclosure. Un-negated, every cross-plate view factor
   came out as exactly 0 (the visibility gate rejected every sample).
2. **Correlated sample positions** (`Radiation.okl`): patch `i`'s and patch
   `j`'s points both derived their stratification cell from the same loop
   index `k`, so they always landed in the "same" relative cell of their
   respective domains — not independent sampling of the double integral.
   This produced a persistent ~1-2% overestimate that did *not* shrink with
   more samples (a converged-but-biased estimator, which is what made it
   findable — statistical noise would have looked different). Fixed by
   hashing `k` to give patch `j` an independent cell.
3. **Three OCCA/OKL parser rejections** (`Radiation.okl`): OCCA's OKL-to-C++
   transpiler (not the downstream compiler) failed to parse a cast of a
   negative literal (`(dfloat)-1e-12`), a combined increment-inside-array-
   index-with-assignment (`stack[sp++] = x`), and a parenthesized cast as a
   bare function-call argument. All are valid C++; rewritten into
   equivalent forms the transpiler accepts (negate-after-cast, split
   increment into its own statement, cast into a named intermediate
   variable before the call).
4. **`nrs->fluid->mesh` → `nrs->meshV`** (`Radiation.cpp`): the former is
   null whenever `[FLUID VELOCITY]` isn't configured (`FLUID=FALSE`),
   which the validation case relies on since it has no flow solver at all.
