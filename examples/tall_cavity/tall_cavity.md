# tall_cavity — Case Files and Problems Solved

A combined natural-convection + surface-radiation example: a tall
(aspect ratio 4) rectangular-prism cavity, one hot wall, one cold wall
directly opposite it, the other four walls held at the hot/cold mean —
every wall no-slip. Exercises the `Radiation::step()` gray-diffuse
radiosity coupling (see `rvf-dev.md`) on a real geometry, alongside a
genuine Navier-Stokes + energy-equation solve.

---

## New files

### `tall_cavity.box`
`genbox` input: `0.08 x 0.08 x 0.32 m` box (4×4×16 elements,
polynomialOrder 7 in `.par`), x0/x1 the hot/cold walls, y0/y1/z0/z1 the
four side walls. Boundary labels are **numeric strings** (`1  `,`2  `,
`3  `) rather than the usual letter codes — required by a build/runtime
problem documented below, not a genbox limitation on its own.

### `tall_cavity.re2`
Generated via `genbox` (see the "getting the mesh built" note below —
genbox's actual output is always literally named `box.re2` and has to
be renamed).

### `tall_cavity.usr`
Legacy Nek5000 Fortran hooks. Every subroutine is the standard empty
stub (copied from `3rd_party/nek5000/core/zero.usr`) **except**
`usrdat2()`, which assigns `boundaryID(ifc,iel)` by matching
`cbc(ifc,iel,1)` against `'1  '`/`'2  '`/`'3  '`. This is the actual
fix for the case's core blocker — see below.

### `tall_cavity.par`
`[GENERAL]` (256-element mesh, `endTime=30`, auto-CFL timestep),
`[FLUID VELOCITY]` (no-slip everywhere, real air properties: ρ=1.2
kg/m³, μ=1.8e-5 Pa·s), `[SCALAR TEMPERATURE]` (all three walls
Dirichlet-fixed via `udfDirichlet`, real air conductivity/ρCp), and
`[RADIATION]` (all 6 walls radiating, `emissivity = 0.9,0.9,0.8` for
hot/cold/sides, `updateFrequency=500`). Uses **dimensional** (Kelvin,
real SI properties) physics throughout, not the nondimensional
Boussinesq convention `rbc.par` uses — deliberately, because a
nondimensional `T∈[0,1]` would be physically meaningless fed into a
real Stefan-Boltzmann flux calculation. The cavity is sized
(8 cm × 8 cm × 32 cm, ΔT=20 K) to land at a computationally tractable
Rayleigh number (Ra≈1.04e6) for this mesh resolution.

### `tall_cavity.udf`
`Radiation::setup()` in `UDF_Setup()`; a Boussinesq buoyancy source
(`userf`, matching `rbc.udf`'s pattern but with real dimensional
`g·β` instead of the pre-baked Ra/Pr trick) added as `nrs->userSource`;
a linear conductive initial temperature profile between the hot and
cold walls; `Radiation::step()` in `UDF_ExecuteStep()`. The
`#ifdef __okl__` block contains `udfDirichlet`, which sets the fixed
wall temperature by `bc->id` (1→310 K, 2→290 K, else→300 K).

Because every wall's temperature is Dirichlet-fixed, the radiative
flux `Radiation::step()` computes has no feedback path into the
temperature solve here — it's a periodic diagnostic of the radiative
exchange between the (convection-independent) fixed wall temperatures,
not a boundary condition. `[RADIATION] updateFrequency=500` is
deliberately coarse since the result never changes between calls.

---

## Problems solved

These were all found only by actually building and running the case
end-to-end — nothing here was visible from reading the code or `.par`
files alone.

### 1. Mixed C++/Fortran toolchain broke the nekRS build
`mpicxx`/`mpicc` on this machine wrap Apple's `clang++` (libc++), while
`mpif77` wraps MacPorts' `gfortran-mp-14` (GNU runtime) — linking
`libnekrs.so` failed with `library 'gfortran' not found`. Pointing the
linker at MacPorts' lib dir directly (`LIBRARY_PATH=...`) "fixed" that
error but broke everything else with hundreds of undefined C++ stdlib
symbols, because it then pulled in MacPorts' `libstdc++` for symbols
that were supposed to resolve via libc++.

**Fix**: reconfigure with a *fully consistent* GCC 14 toolchain for
C/C++ *and* Fortran, via `MPICH_CC`/`MPICH_CXX`/`MPICH_FC` env vars
(which redirect the `mpicc`/`mpicxx`/`mpif77` wrappers to
`gcc-mp-14`/`g++-mp-14`/`gfortran-mp-14`), then a full clean rebuild.
These same three env vars are also required *at case-run time* (not
just at build time) — nekRS compiles a small nek5000 C interface layer
on first run of each case, using `CC` from the installed
`nekrs.conf`, which again resolves through the `mpicc` wrapper.

### 2. GCC + macOS SDK headers: `_Alignof` incompatibility
GCC on this macOS SDK chokes on `_Alignof` used in system headers
(`mach/arm/_structs.h`) — `error: '_Alignof' was not declared in this
scope`. Hit this twice, in two different places:
- Building nekRS itself (`std::filesystem` detection, and generally
  anything including `<cstdlib>`/similar). **Fix**: `CXXFLAGS=
  -D_Alignof=alignof` — **C++ only**. Setting the same on `CFLAGS`
  breaks plain-C third-party libraries (HYPRE, ADIOS2's C
  dependencies) instead, since C doesn't have `alignof` as a keyword.
- OCCA's *runtime* kernel JIT compiler, which is configured
  independently via nekRS's own `OCCA_CXXFLAGS` CMake variable (not
  inherited from the main build's `CXXFLAGS`) and gets baked into the
  installed `nekrs.conf`. **Fix**: also pass
  `-DOCCA_CXXFLAGS="... -D_Alignof=alignof"` at configure time.

### 3. genbox never writes numeric boundary IDs for ordinary faces (the real blocker)
The mesh loaded with `NboundaryIDs: 0`, which then aborted at
`fluidSolver.cpp:110` ("`boundaryTypeMap` (3) does not match number of
boundary IDs in mesh (0)"). Traced (with help from an
investigation agent that read the actual Fortran source and confirmed
it against a byte-level parse of a shipped example's `.re2`) to:

- genbox writes the 3-character CBC label into the `.re2`, but leaves
  the 5 real-valued "BC parameter" slots at zero for any label other
  than periodic (`'P  '`) — regardless of whether the label is a
  letter code or a numeric string. Confirmed by decoding a boundary
  record from `examples/channel/channel.re2`: `bc(1..5) = [0,0,0,0,0]`
  for a `'SYM'` face.
- nekRS's automatic CBC-label→numeric-ID import (`gen_bcmap()`, which
  *would* have derived IDs 1/2/3 from the labels automatically) only
  runs when `useNek()` is true — and `useNek()` is set permanently
  false, for *every* field, the moment *any* field gets an explicit
  `boundaryTypeMap` in `.par` (`bdryBase::setupField`,
  `bdryBase.hpp:130-137`). This case needs `boundaryTypeMap` on both
  `[FLUID VELOCITY]` and `[SCALAR TEMPERATURE]` (different BC types
  needed per field on the same raw IDs), so the automatic path is
  unavailable.
- `[MESH] boundaryIDMap` (documented in `doc/parHelp.txt`) does **not**
  fix this — it only *remaps* already-positive raw IDs; there's
  nothing positive to remap.
- The actual, only-documented-by-example fix: three shipped examples
  (`examples/channel/channel.usr`, `examples/mv_cyl/mv_cyl.usr`,
  `examples/periodicHill/periodicHill.usr`) all assign `boundaryID`
  by hand in a `usrdat2()` Fortran subroutine, matching against the
  CBC text (which *does* survive the round trip through the `.re2`
  file, unlike the numeric slot). `tall_cavity` had no `.usr` file at
  all — that was the actual bug. Added one (see above).

### 4. `udfDirichlet`/`udfNeumann` must live in `.udf`, not `.oudf`
A separate `.oudf` file is for extra standalone `@kernel`s (like
`gabls1`'s `coriolis` kernel) — nekRS's UDF-completeness check
(`adjustOudf` in `src/core/udf/udf.cpp`) only scans the preprocessed
`#ifdef __okl__` block *inside* the `.udf` file for `udfDirichlet`/
`udfNeumann`/`udfRobin`. An initial `tall_cavity.oudf` with
`udfDirichlet` in it was silently never scanned, producing "Cannot
find required okl function udfDirichlet!" even though the function
existed, just in the wrong file. Fixed by moving it into
`tall_cavity.udf`'s `#ifdef __okl__` block (matching `rbc.udf`'s
pattern) and deleting the `.oudf`.

---

## Getting the mesh built

`genbox` always writes its output to a file literally named `box.re2`
regardless of the `.box` file's declared name or the name typed at its
interactive prompt — rename it afterward:

```
echo "tall_cavity.box" | genbox
mv box.re2 tall_cavity.re2
```

## Result

The case runs to completion (`t=0→30 s`, ~3000 steps), stable
throughout (CFL bounded ~1.1–2.2, all residuals converging every
step, no divergence). `Radiation::setup()` finds 288 radiating patches
(64 each for hot/cold, 64+64 for the two larger side walls, 16+16 for
the two smaller ones — matches the hand-computed element-face count
exactly). View factors from `tall_cavity_radiation_viewfactors_groups.csv`:
`F_12 = F_21 = 0.346` (hot↔cold), `F_13 = 0.657`, `F_23 = 0.651`
(reciprocity `A_i F_ij = A_j F_ji` holds exactly; each group's row
sums to ~1.003, i.e. ~0.3% Monte Carlo noise at `nSamples=8192`).

The periodic radiative-flux diagnostic (identical at every report,
as expected since all wall temperatures are fixed) settles to:

| bID | wall  | T (K) | net radiative flux (W/m², avg over patches) |
|-----|-------|-------|-----------------------------------------------|
| 1   | hot   | 310   | +69.1 (net radiator)                          |
| 2   | cold  | 290   | −66.6 (net absorber)                          |
| 3   | sides | 300   | −5.1 (net absorber)                           |

The side walls being a net *absorber* despite sitting exactly at the
hot/cold arithmetic mean is a good physical sanity check, not an
anomaly: blackbody emissive power scales as T⁴ (convex), so
`E_b(310) + E_b(290) > 2·E_b(300)` — a mean-temperature surface
necessarily receives more radiation from a hot/cold pair than it would
emit at that same mean temperature.
