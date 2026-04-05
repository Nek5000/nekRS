# nekRS Boundary Condition Reference

## bcData Struct Fields

Available inside `udfDirichlet`, `udfNeumann`, `udfRobin` functions:

### Position and Time
```
bc->x, bc->y, bc->z     // Boundary node coordinates
bc->time                 // Current simulation time
```

### Surface Geometry
```
bc->nx, bc->ny, bc->nz   // Unit outward normal
bc->t1x, bc->t1y, bc->t1z // Tangent vector 1
bc->t2x, bc->t2y, bc->t2z // Tangent vector 2
```

### Fluid Fields (set these in udfDirichlet)
```
bc->uxFluid, bc->uyFluid, bc->uzFluid  // Velocity Dirichlet values
bc->pFluid                              // Pressure value
bc->uxFluidInt, bc->uyFluidInt, bc->uzFluidInt  // Interior interpolated values
```

### Scalar Fields
```
bc->sScalar      // Scalar Dirichlet value (set in udfDirichlet)
bc->fluxScalar   // Scalar normal flux (set in udfNeumann)
bc->sScalarInt   // Interior interpolated scalar
bc->sInfScalar   // Far-field scalar (Robin BC)
bc->h            // Heat transfer coefficient (Robin BC)
bc->diffCoeff    // Diffusion coefficient at boundary
bc->transCoeff   // Transport coefficient at boundary
bc->idScalar     // Scalar index
```

### Traction (Neumann for velocity)
```
bc->tr1, bc->tr2   // Tangential traction components
```

### Mesh Motion (Geometry)
```
bc->uxGeom, bc->uyGeom, bc->uzGeom  // Mesh velocity Dirichlet values
```

### Metadata
```
bc->id           // Boundary tag ID (from mesh)
bc->idxVol       // Volume node index
bc->fieldOffset  // Field array offset
bc->usrwrk       // Read-only pointer to user workspace array
```

### Field Identification
```cpp
isField("fluid velocity")      // Check if current field is velocity
isField("fluid pressure")      // Check if current field is pressure
isField("scalar temperature")  // Named scalar
isField("scalar00")            // Scalar by index
isField("geom")                // Mesh motion field
```

---

## Function Templates

### udfDirichlet — Set Dirichlet Values

```okl
void udfDirichlet(bcData *bc)
{
  // Velocity inlet: parabolic profile
  if (isField("fluid velocity")) {
    if (bc->id == 1) {  // inlet boundary
      const dfloat y = bc->y;
      const dfloat H = 1.0;
      bc->uxFluid = 4.0 * y * (H - y) / (H * H);
      bc->uyFluid = 0.0;
      bc->uzFluid = 0.0;
    }
  }

  // Temperature: fixed wall temperature
  if (isField("scalar temperature")) {
    if (bc->id == 1) {
      bc->sScalar = 1.0;  // hot wall
    } else if (bc->id == 2) {
      bc->sScalar = 0.0;  // cold wall
    }
  }
}
```

### udfNeumann — Set Normal Flux

```okl
void udfNeumann(bcData *bc)
{
  if (isField("scalar temperature")) {
    // Constant heat flux
    bc->fluxScalar = 1.0;
  }
}
```

### udfRobin — Set Convective/Mixed BC

```okl
void udfRobin(bcData *bc)
{
  if (isField("scalar temperature")) {
    // h*(T - T_inf) + k*dT/dn = 0
    bc->h = 10.0;            // heat transfer coefficient
    bc->sInfScalar = 300.0;  // ambient temperature
    // bc->diffCoeff is set automatically from .par
  }
}
```

---

## Common BC Recipes

### Parabolic Inlet (2D channel)
```okl
if (isField("fluid velocity")) {
  const dfloat y = bc->y;
  bc->uxFluid = 6.0 * y * (1.0 - y);  // normalized for unit height
  bc->uyFluid = 0.0;
  bc->uzFluid = 0.0;
}
```

### Poiseuille Pipe Inlet
```okl
if (isField("fluid velocity")) {
  const dfloat r2 = bc->y * bc->y + bc->z * bc->z;
  const dfloat R = 0.5;
  bc->uxFluid = 2.0 * U_bulk * (1.0 - r2 / (R * R));
  bc->uyFluid = 0.0;
  bc->uzFluid = 0.0;
}
```

### Time-Varying Inlet
```okl
if (isField("fluid velocity")) {
  const dfloat amp = 0.1;
  const dfloat freq = 2.0;
  bc->uxFluid = 1.0 + amp * sin(freq * bc->time);
  bc->uyFluid = 0.0;
  bc->uzFluid = 0.0;
}
```

### Temperature Ramp
```okl
if (isField("scalar temperature")) {
  const dfloat t_ramp = 1.0;  // ramp duration
  const dfloat T_hot = 1.0;
  const dfloat factor = tanh(bc->time / t_ramp);
  bc->sScalar = T_hot * factor;
}
```

### Rotating Wall
```okl
if (isField("fluid velocity")) {
  const dfloat omega = 1.0;
  bc->uxFluid = -omega * bc->y;
  bc->uyFluid =  omega * bc->x;
  bc->uzFluid = 0.0;
}
```

### Using usrwrk for Computed BCs
```okl
// In .oudf - read data set from UDF side
void udfDirichlet(bcData *bc)
{
  if (isField("fluid velocity")) {
    const int idxU = 0;  // usrwrk slot indices
    const int idxV = 1;
    const int idxW = 2;
    bc->uxFluid = bc->usrwrk[bc->idxVol + idxU * bc->fieldOffset];
    bc->uyFluid = bc->usrwrk[bc->idxVol + idxV * bc->fieldOffset];
    bc->uzFluid = bc->usrwrk[bc->idxVol + idxW * bc->fieldOffset];
  }
}
```

---

## Boundary Type Mapping Rules

1. Types in `boundaryTypeMap` are **comma-separated** and **positional**: first = boundary ID 1, second = boundary ID 2, etc.
2. The number of entries must match the number of unique boundary IDs in the `.re2` mesh.
3. Compound types use `/`: `zeroDirichletN/zeroNeumann` means zero normal velocity + zero tangential gradient.
4. Every `udfDirichlet` / `udfNeumann` / `udfRobin` boundary must have a corresponding handler in the `.oudf` or `__okl__` section.
5. Periodic boundaries are handled by the mesh (not in `boundaryTypeMap`) — they don't appear in the map.
