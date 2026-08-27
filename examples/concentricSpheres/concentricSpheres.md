# concentricSpheres — Case Files, Pipeline, and Results

A pure-radiation diagnostic case: two concentric spheres (inner
r=0.5 m at 1000 K, outer r=1.0 m at 500 K, both ε=0.85), reproducing
Section 4.2 ("Concentric spheres") of Yuan, Dai, Shaver, Gottems, and
Merzari, *"Implementation of the Surface-to-Surface Thermal Radiation
Solver in Spectral Element CFD code NekRS"*.

Unlike `tall_cavity`/`radiationPlates`, `genbox` can't produce curved
geometry, so the mesh was built with `gmsh` (built-in kernel, `.geo`
script) and converted with `gmsh2nek`, following the same all-Dirichlet
radiation-diagnostic pattern otherwise (no flow of interest, fixed
wall temperatures, `Radiation::step()` called every step).

---

## New files

### `make_mesh.py`
Generates `concentricSpheres.geo`: a **cubed-sphere** hex mesh of the
shell. genbox can't do curved geometry, and a naive single-block
octant of a shell is topologically a triangular prism (6 vertices, 5
faces) — not a hexahedron nekRS/Nek5000 elements require. The
standard resolution is the cubed sphere: the 6 cube faces each become
one curved hex "panel", bounded by 2 spherical-cap surfaces (built
from exact great-circle arcs, `Circle(id) = {pStart, origin, pEnd};`
— curvature captured exactly, not approximated) and 4 radial side
faces shared with the adjacent panels. All 6 panels reuse the same
corner/edge Point/Curve IDs by construction, so the mesh is exactly
conformal with no separate coincidence-merge step. `Physical Surface`
groups tag the inner cap (`"hot"`, ID 1) and outer cap (`"cold"`, ID 2).

Parameters: `N=4` points per panel edge (3 elements/edge, 9/panel,
54 total), `NR=2` (1 element radially — the radial direction is
straight, no curvature-fidelity concern).

### `concentricSpheres.geo`, `.msh`, `.re2`
Generated in sequence: `python3 make_mesh.py` → `gmsh concentricSpheres.geo
-3 -format msh22 -o concentricSpheres.msh -nopopup` → `gmsh2nek`
(interactive: `3` / `concentricSpheres` / `0` / `0` / `concentricSpheres`
— dimension, fluid mesh name, no solid mesh, no periodic pairs, output
name).

### `concentricSpheres.usr`
Plain `zero.usr` template, **no boundary-ID fix needed**. Unlike
genbox (which always zeros the `.re2`'s numeric BC-parameter slot for
ordinary faces, requiring the `usrdat2()` workaround in `tall_cavity`
and `radiationPlates`), `gmsh2nek` writes the Physical Surface tag
directly into `bc(5,...)` (`gmsh2nek.f90:1551`), so nekRS's boundary-ID
import picks it up with no manual assignment.

### `concentricSpheres.par` / `.udf`
Same all-Dirichlet diagnostic structure as `radiationPlates`:
`[FLUID VELOCITY]` no-slip (uninteresting, nothing forces it),
`[SCALAR TEMPERATURE]` `udfDirichlet` (1→1000 K, 2→500 K),
`[RADIATION]` both surfaces radiating **and** both listed as
`obstructionBoundaryIDs` (see Problem 4 below).

---

## Problems solved

### 1. A leftover gmsh GUI option silently multiplied every element by 8×
First mesh attempt (before any of the geometry was even suspect) gave
48 hex elements instead of the intended 6. Traced to
`~/.gmsh-options` — a **persisted setting from an unrelated prior gmsh
GUI session** (`Mesh.SubdivisionAlgorithm = 2`, "all hexahedra"),
auto-loaded on every gmsh invocation per that file's own header
comment, silently overriding the `.geo` script's transfinite
subdivision with a global 2×2×2 refinement. Confirmed via a minimal
1-hex reproduction case (8 hexes instead of 1) and fixed by adding
`Mesh.SubdivisionAlgorithm = 0;` explicitly inside the generated
`.geo`, rather than relying on ambient user state.

### 2. `gmsh2nek` needs a `$PhysicalNames` section — numeric-only tags don't produce one
`gmsh2nek.f90:871-892` loops reading lines until it finds
`$PhysicalNames`, and crashes with a Fortran "End of file" error if
that section is absent. `Physical Surface(1) = {...}` (tag number
only, no name string) does not make gmsh write that section; using
named groups (`Physical Surface("hot", 1) = {...};`) does. Confirmed
by inspecting the raw `.msh` output before and after the fix.

### 3. `gmsh2nek`'s interactive prompt sequence (undocumented in the README)
Traced from `gmsh2nek.f90`'s main program: `Enter mesh dimension` (3)
→ (inside `read_input_name`, called immediately after) `Input fluid
.msh file name` (base name, no extension — the code appends `.msh`
itself) → `Do you have solid mesh?` (0) → `Enter number of periodic
boundary surface pairs` (0) → `please give re2 file name` (base name
again). Getting this sequence wrong (e.g. skipping the periodicity
prompt) corrupts the whole read and produces a "Bad integer for item
1 in list input" crash several prompts later, not at the point of the
actual mistake.

### 4. Coarse single-element-per-panel geometry: sphere surface area 12% low
The first working mesh (`N=2`, 1 huge ~90°-wide quadratic element per
panel) ran cleanly with no errors, but the reported patch areas were
**12% below the true sphere area** (`A_1=2.767` vs. true `π=3.1416`
for `r=0.5`; same ratio for the outer sphere). A 20-node quadratic
element only controls curvature along its edges, not its face
interior — filling a ~90° curved patch from just boundary information
under-represents the true bulge. Fixed by refining to `N=4` (3×3
elements per panel face, 30° each): area error dropped to 0.22%,
consistent with the expected roughly-`θ⁴` scaling of quadratic-element
curvature error.

### 5. `obstructionBoundaryIDs` is not automatic — radiating surfaces don't occlude each other by default
First full run with the refined mesh gave `F_21=0.2499`,
`F_22=0.9999` for the outer (concave) sphere — a **row sum of 1.25**,
which is impossible for genuine view factors (energy leaving a
surface can't sum to more than 1). Traced (with help from a targeted
code-reading pass over `Radiation.cpp`/`Radiation.okl`) to: the
occlusion BVH is only built when `obstructionBoundaryIDs` is
non-empty (`Radiation.cpp`: `hasObstruction = !obstructionIDs.empty()`),
and the OKL kernel skips the occlusion test entirely when no BVH
exists (`Radiation.okl`: `blocked=0` unless `nBvhNodes>0`). A surface
listed only in `radiatingBoundaryIDs` does **not** automatically also
act as an occluder for other radiating surfaces' mutual visibility —
that's a real, load-bearing distinction, not a bug, but it's easy to
miss for a concave-surface case like this one (the outer sphere's
self-view genuinely needs to be blocked by the inner sphere). Fixed
by adding `obstructionBoundaryIDs = 1,2` to `[RADIATION]`. After the
fix: row sum for the outer sphere = 0.9971 (0.3% MC noise, matching
the pattern already established in `tall_cavity`/`radiationPlates`).

---

## Results

Measured view factors (`concentricSpheres_radiation_viewfactors_groups.csv`,
after fix #4 and #5):

```
F_11 = 0            (inner sphere is convex, can't see itself -- exact)
F_12 = 0.999121      (inner -> outer)
F_21 = 0.249782      (outer -> inner; exact analytic value is A_1/A_2 = 0.25 -- 0.09% off)
F_22 = 0.747268      (outer self-view, occluded by the inner sphere)
```

Reciprocity (`A_1 F_12 = A_2 F_21`): `3.1347×0.99912 = 3.1319` vs.
`12.5386×0.24978 = 3.1318` — matches to 0.005%. Row sums: 0.9991
(inner) and 0.9971 (outer), both within ~0.3% MC noise
(`nSamples=16384`).

Radiative flux, compared to the paper's Eq. 12 analytical formula
(evaluated at our own measured/exact geometry, `A_1/A_2 = 0.25`):

| surface | T (K) | nekRS (this case) | Eq. 12 analytical | paper Table II (different geometry) |
|---------|-------|--------------------|--------------------|----------------------------------------|
| inner (hot)  | 1000 | +43,560.4 W/m² | +43,552.6 W/m² | +4.3522×10⁴ (analytical 4.36×10⁴) |
| outer (cold) | 500  | −10,861.5 W/m² | — | −1.105×10⁴ (analytical −1.09×10⁴) |

**nekRS vs. its own case's analytical formula: 0.018% agreement** —
essentially exact, and tighter than the parallel-plates case, because
concentric spheres is precisely the geometry Eq. 12/13 assumes (no
aspect-ratio approximation needed, unlike `radiationPlates`).

Two independent checks:
- **Global energy balance**: `A_1 q_1 + A_2 q_2 = 359.7`, which is
  0.26% of `A_1 q_1` — closes to within MC noise.
- **Lumped 2-node network** (solving Eqs. 1–4 by hand with the
  measured `F` matrix) reproduces nekRS's own patch-resolved output
  almost exactly (43,560.4 vs. 43,560.4; −10,863.3 vs. −10,861.5).
  Unlike `radiationPlates`'s side walls, every patch within a group
  here is geometrically equivalent by spherical symmetry, so the
  lumped and patch-resolved solutions agree — a good internal
  consistency check in its own right.
