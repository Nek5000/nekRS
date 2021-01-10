#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <gencon-impl.h>
#include <parRSB.h>

void fparRSB_findConnectivity(long long *vertexId, double *coord, int *nelt, int *ndim,
                              long long *periodicInfo, int *nPeriodicFaces, double *tol, MPI_Fint *fcomm,
                              int *verbose, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*fcomm);
  *err = parRSB_findConnectivity(vertexId, coord, *nelt, *ndim, periodicInfo, *nPeriodicFaces,
                                 *tol, c, *verbose);
}

// coord[nelt,nv,ndim] - in, vertices are orders in preprocessor ordering
// vertexid[nelt,nv] - out
int parRSB_findConnectivity(long long *vertexid, double *coord, int nelt, int ndim, long long *periodicInfo,
                            int nPeriodicFaces, double tol, MPI_Comm comm, int verbose)
{
  struct comm c;
  comm_init(&c, comm);
  uint rank = c.id;
  uint size = c.np;

  if (rank == 0 && verbose > 0)
    printf("Generating connectivity ...");
  fflush(stdout);

  if (verbose > 0)
    printf("\n\tnelt/ndim/nperiodic: %d/%d/%d\n", nelt, ndim, nPeriodicFaces);

  double t_con = 0.0;
  comm_barrier(&c);
  t_con -= comm_time();

  Mesh mesh;
  mesh_init(&mesh, nelt, ndim);

  slong out[2][1], buff[2][1];
  slong in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, buff);
  ulong start = out[0][0];
  ulong nelgt = mesh->nelgt = mesh->nelgv = out[1][0];

  int nelt_ = nelgt/size;
  int nrem = nelgt - nelt_*size;
  if (rank >= (size - nrem))
    nelt_++;
  assert(nelt==nelt_);

  int nvertex = mesh->nVertex;
  uint nunits = nvertex*nelt;

  struct Point_private p;
  uint i, j, k, l;
  for (i = 0; i < nelt; i++){
    for (k = 0; k < nvertex; k++){
      j = PRE_TO_SYM_VERTEX[k];
      for(l = 0; l < ndim; l++)
        p.x[l] = coord[i*nvertex*ndim + j*ndim + l];
      p.elementId = start + i;
      p.sequenceId = nvertex*(start + i) + k;
      p.origin = rank;

      array_cat(struct Point_private, &mesh->elements, &p, 1);
    }
  }
  assert(mesh->elements.n == nunits);

  struct Boundary_private b;
  for (i = 0; i<nPeriodicFaces; i++) {
    b.elementId = periodicInfo[4*i+0]-1;
    b.faceId = PRE_TO_SYM_FACE[periodicInfo[4*i+1]-1];
    b.bc[0] = periodicInfo[4*i+2]-1;
    b.bc[1] = PRE_TO_SYM_FACE[periodicInfo[4*i+3]-1];
    array_cat(struct Boundary_private,&mesh->boundary,&b,1);
  }
  assert(mesh->boundary.n==nPeriodicFaces);

  transferBoundaryFaces(mesh,&c);

  findMinNeighborDistance(mesh);
  findSegments(mesh,&c,tol,1);
  setGlobalID(mesh,&c);
  sendBack(mesh,&c);
  faceCheck(mesh,&c);
  matchPeriodicFaces(mesh,&c);

  // copy output
  Point ptr = mesh->elements.ptr;
  for (i = 0; i < nelt*nvertex; i++)
    vertexid[i] = ptr[i].globalId + 1;

  comm_barrier(&c);
  t_con += comm_time();

  if (rank == 0 && verbose > 0)
    printf("\tfinished in %g s\n",t_con);

  mesh_free(mesh);
  comm_free(&c);

  return 0;
}
