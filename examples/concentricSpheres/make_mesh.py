#!/usr/bin/env python3
"""
Generates a "cubed sphere" hex mesh of the shell between two concentric
spheres (r_in, r_out), for the Yuan et al. Sec 4.2 concentric-sphere
radiation benchmark. genbox can't do curved geometry, so this writes a
gmsh .geo script (built-in kernel, not OCC): each of the 6 cube faces
becomes one curved hex block ("panel"), bounded by two spherical caps
(built from exact great-circle arcs, so the curvature is captured
exactly, not approximated) and 4 radial side faces shared with the
adjacent panels. All 6 panels share their corner/edge geometry by
construction (same Point/Line IDs reused, not just numerically close),
so the mesh is exactly conformal -- no separate merge/tolerance step
needed.

Physical Surface 1 = all 6 inner-sphere caps (hot), 2 = all 6
outer-sphere caps (cold). gmsh2nek writes these directly into the
.re2's boundaryID (bc(5,...)) -- unlike genbox, which never writes a
usable numeric ID for ordinary faces, gmsh-tagged meshes don't need
the usrdat2() boundaryID workaround examples/tall_cavity and
examples/radiationPlates both needed.
"""
import itertools

R_IN = 0.5
R_OUT = 1.0
N = 4  # points per panel edge (>=2); N=4 -> 3 elements per panel edge (9/panel, 54 total per radial layer)
NR = 2  # points radially (>=2); NR=2 -> 1 element radially
# N=2 (1 huge ~90-deg element per panel) was tried first and measured
# ~12% low on sphere surface area (2.767 vs true pi=3.1416 for r=0.5)
# -- a single quadratic (HEXA20) element can't represent that much
# curvature; refining to N=4 (30-deg elements) is expected to cut the
# area error by ~O(theta^4), i.e. roughly 80x, per patch.

lc = 1.0  # unused (transfinite everywhere), required by Point syntax

lines = []
def emit(s):
    lines.append(s)

emit(f"// Cubed-sphere shell mesh: r_in={R_IN}, r_out={R_OUT}")
emit("SetFactory(\"Built-in\");")
emit("")

# --- Origin ---
emit(f"Point(1) = {{0, 0, 0, {lc}}};")
ORIGIN = 1
next_pt = 2

# --- 8 corner directions, 2 radii each ---
signs = list(itertools.product([-1, 1], repeat=3))  # (sx,sy,sz)
corner_pt = {}  # (sx,sy,sz, 'in'|'out') -> point id
import math
sqrt3 = math.sqrt(3)
for (sx, sy, sz) in signs:
    ux, uy, uz = sx / sqrt3, sy / sqrt3, sz / sqrt3
    for radius, tag in ((R_IN, 'in'), (R_OUT, 'out')):
        pid = next_pt
        next_pt += 1
        corner_pt[(sx, sy, sz, tag)] = pid
        emit(f"Point({pid}) = {{{ux*radius:.10f}, {uy*radius:.10f}, {uz*radius:.10f}, {lc}}};")
emit("")

# --- 12 cube edges: fixed axis pair, varying axis ---
# an edge is identified by which two axes are FIXED (with signs) and which axis varies
def other_axes(axis):
    return [a for a in range(3) if a != axis]

edge_arc = {}  # (frozenset of 2 fixed (axis,sign) pairs) -> {'in':lineid,'out':lineid}
next_line = 1

def corner_from_fixed(varying_axis, varying_sign, fixed):
    # fixed: dict axis->sign for the two fixed axes
    s = [0, 0, 0]
    s[varying_axis] = varying_sign
    for a, sgn in fixed.items():
        s[a] = sgn
    return tuple(s)

edges = []  # list of (varying_axis, fixed_axis0, fixed_sign0, fixed_axis1, fixed_sign1)
for varying_axis in range(3):
    fa0, fa1 = other_axes(varying_axis)
    for s0, s1 in itertools.product([-1, 1], repeat=2):
        edges.append((varying_axis, fa0, s0, fa1, s1))

edge_key_to_arc = {}
for (va, fa0, s0, fa1, s1) in edges:
    fixed = {fa0: s0, fa1: s1}
    cA = corner_from_fixed(va, -1, fixed)
    cB = corner_from_fixed(va, +1, fixed)
    key = (va, fa0, s0, fa1, s1)
    arcs = {}
    for tag in ('in', 'out'):
        pA = corner_pt[cA + (tag,)]
        pB = corner_pt[cB + (tag,)]
        lid = next_line
        next_line += 1
        emit(f"Circle({lid}) = {{{pA}, {ORIGIN}, {pB}}};")
        arcs[tag] = lid
    edge_key_to_arc[key] = (cA, cB, arcs)
emit("")

# --- 8 radial lines (corner in -> corner out) ---
radial_line = {}
for s in signs:
    pin = corner_pt[s + ('in',)]
    pout = corner_pt[s + ('out',)]
    lid = next_line
    next_line += 1
    emit(f"Line({lid}) = {{{pin}, {pout}}};")
    radial_line[s] = lid
emit("")

# --- Transfinite on all curves ---
for key, (cA, cB, arcs) in edge_key_to_arc.items():
    emit(f"Transfinite Curve {{{arcs['in']}}} = {N};")
    emit(f"Transfinite Curve {{{arcs['out']}}} = {N};")
for s, lid in radial_line.items():
    emit(f"Transfinite Curve {{{lid}}} = {NR};")
emit("")

# --- 6 panels ---
# panel defined by fixed axis + sign; corners ordered CCW as seen from outside
panel_defs = []
for axis in range(3):
    for sign in (-1, 1):
        panel_defs.append((axis, sign))

def panel_corners_ccw(axis, sign):
    # the two free axes, in an order chosen so traversal is a consistent CCW loop
    free = other_axes(axis)
    a0, a1 = free
    # order: (-,-),(+,-),(+,+),(-,+) in (a0,a1) - consistent quad loop
    order = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    corners = []
    for (s0, s1) in order:
        s = [0, 0, 0]
        s[axis] = sign
        s[a0] = s0
        s[a1] = s1
        corners.append(tuple(s))
    return corners, a0, a1

def find_edge_key(cA, cB):
    # cA,cB differ in exactly one coordinate; find matching edge key (order-independent)
    diffs = [i for i in range(3) if cA[i] != cB[i]]
    assert len(diffs) == 1
    va = diffs[0]
    fa0, fa1 = other_axes(va)
    fixed = {fa0: cA[fa0], fa1: cA[fa1]}
    key = (va, fa0, fixed[fa0], fa1, fixed[fa1])
    return key

surf_id = 1
cap_surfaces = {'in': [], 'out': []}
radial_surf_cache = {}  # edge key -> surface id (shared between 2 panels)

surface_lines = []

for (axis, sign) in panel_defs:
    corners, a0, a1 = panel_corners_ccw(axis, sign)
    # --- caps (r_in, r_out) ---
    for tag in ('in', 'out'):
        loop_lines = []
        for i in range(4):
            cA = corners[i]
            cB = corners[(i + 1) % 4]
            key = find_edge_key(cA, cB)
            cAk, cBk, arcs = edge_key_to_arc[key]
            lid = arcs[tag]
            # orientation: if our traversal cA->cB matches the arc's own cAk->cBk, positive; else negate
            sign_l = lid if cA == cAk else -lid
            loop_lines.append(sign_l)
        sid = surf_id
        surf_id += 1
        cl_id = 10000 + sid
        surface_lines.append(f"Curve Loop({cl_id}) = {{{','.join(str(x) for x in loop_lines)}}};")
        surface_lines.append(f"Surface({sid}) = {{{cl_id}}};")
        pts_in_order = [corner_pt[c + (tag,)] for c in corners]
        surface_lines.append(f"Transfinite Surface {{{sid}}} = {{{','.join(str(p) for p in pts_in_order)}}};")
        surface_lines.append(f"Recombine Surface {{{sid}}};")
        cap_surfaces[tag].append(sid)

    # --- 4 radial side surfaces ---
    for i in range(4):
        cA = corners[i]
        cB = corners[(i + 1) % 4]
        key = find_edge_key(cA, cB)
        if key in radial_surf_cache:
            continue  # already built by the neighboring panel
        cAk, cBk, arcs = edge_key_to_arc[key]
        arc_in = arcs['in']
        arc_out = arcs['out']
        lineA = radial_line[cAk]  # cAk_in -> cAk_out
        lineB = radial_line[cBk]
        # loop: arc_in(cAk->cBk), lineB(cBk_in->cBk_out), -arc_out(cBk->cAk), -lineA(cAk_out->cAk_in)
        sid = surf_id
        surf_id += 1
        cl_id = 10000 + sid
        surface_lines.append(f"Curve Loop({cl_id}) = {{{arc_in},{lineB},{-arc_out},{-lineA}}};")
        surface_lines.append(f"Surface({sid}) = {{{cl_id}}};")
        pin_A = corner_pt[cAk + ('in',)]
        pin_B = corner_pt[cBk + ('in',)]
        pout_B = corner_pt[cBk + ('out',)]
        pout_A = corner_pt[cAk + ('out',)]
        surface_lines.append(f"Transfinite Surface {{{sid}}} = {{{pin_A},{pin_B},{pout_B},{pout_A}}};")
        surface_lines.append(f"Recombine Surface {{{sid}}};")
        radial_surf_cache[key] = sid

emit("\n".join(surface_lines))
emit("")

# --- 6 volumes ---
vol_lines = []
for vidx, (axis, sign) in enumerate(panel_defs):
    corners, a0, a1 = panel_corners_ccw(axis, sign)
    # find this panel's cap surface ids (recompute by regenerating in the same call order)
    pass

# Need the cap surface id assigned per panel; redo bookkeeping with a dict this time
emit_recap = []
surf_id = 1
panel_cap_sid = {}
radial_surf_cache2 = {}
for (axis, sign) in panel_defs:
    corners, a0, a1 = panel_corners_ccw(axis, sign)
    for tag in ('in', 'out'):
        panel_cap_sid[(axis, sign, tag)] = surf_id
        surf_id += 1
    for i in range(4):
        cA = corners[i]
        cB = corners[(i + 1) % 4]
        key = find_edge_key(cA, cB)
        if key not in radial_surf_cache2:
            radial_surf_cache2[key] = surf_id
            surf_id += 1

for (axis, sign) in panel_defs:
    corners, a0, a1 = panel_corners_ccw(axis, sign)
    sid_in = panel_cap_sid[(axis, sign, 'in')]
    sid_out = panel_cap_sid[(axis, sign, 'out')]
    side_sids = []
    for i in range(4):
        cA = corners[i]
        cB = corners[(i + 1) % 4]
        key = find_edge_key(cA, cB)
        side_sids.append(radial_surf_cache2[key])
    sl_id = 20000 + sid_in
    vol_lines.append(f"Surface Loop({sl_id}) = {{{sid_in},{sid_out},{','.join(str(s) for s in side_sids)}}};")
    vid = sid_in  # reuse numbering space, unique per panel
    vol_lines.append(f"Volume({vid}) = {{{sl_id}}};")
    pts_in = [corner_pt[c + ('in',)] for c in corners]
    pts_out = [corner_pt[c + ('out',)] for c in corners]
    vol_lines.append(f"Transfinite Volume {{{vid}}} = {{{','.join(str(p) for p in pts_in + pts_out)}}};")
    vol_lines.append(f"Recombine Volume {{{vid}}};")

emit("\n".join(vol_lines))
emit("")

all_vols = [panel_cap_sid[(axis, sign, 'in')] for (axis, sign) in panel_defs]
# Named physical groups: gmsh2nek (gmsh2nek.f90:871-892) requires a
# $PhysicalNames section in the .msh to find boundary tags at all --
# numeric-only "Physical Surface(1) = {...}" (no name string) does NOT
# get gmsh to emit that section, and gmsh2nek then hits EOF searching
# for it. Discovered empirically: the unnamed form produced a .msh
# gmsh2nek couldn't read.
emit(f"Physical Surface(\"hot\", 1) = {{{','.join(str(panel_cap_sid[(a,s,'in')]) for (a,s) in panel_defs)}}};")
emit(f"Physical Surface(\"cold\", 2) = {{{','.join(str(panel_cap_sid[(a,s,'out')]) for (a,s) in panel_defs)}}};")
emit(f"Physical Volume(\"fluid\", 1) = {{{','.join(str(v) for v in all_vols)}}};")
emit("")
emit("Mesh.ElementOrder = 2;")
emit("Mesh.SecondOrderIncomplete = 1;")
emit("Mesh.RecombineAll = 1;")
# Override any ambient ~/.gmsh-options state (e.g. a persisted
# Mesh.SubdivisionAlgorithm from an unrelated GUI session) that would
# silently subdivide every transfinite element -- observed empirically:
# without this, Gmsh 4.14.0 on this machine turned each single-hex
# panel into 8 via a leftover "all hexahedra" subdivision setting.
emit("Mesh.SubdivisionAlgorithm = 0;")
emit("Coherence;")

with open("concentricSpheres.geo", "w") as f:
    f.write("\n".join(lines) + "\n")

print("wrote concentricSpheres.geo")
print(f"points: {next_pt-1}, curves: {next_line-1}, surfaces: {surf_id-1}, volumes: {len(all_vols)}")
