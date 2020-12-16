/*
Parition mesh using Nek5000's vertex connectivity (con) file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <gencon.h>
#include <genmap.h>
#include <parRSB.h>

#include "quality.h"

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {
  int proc;
  long long vtx[MAXNV];
} elm_data;

#define EXIT_ERROR() do{\
  MPI_Finalize();\
  return EXIT_FAILURE;\
} while(0)

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  if(argc!=3){
    if(rank==0)
      printf("Usage: %s <#nread> <mesh file>\n",argv[0]);
    EXIT_ERROR();
  }

  int n_read = atoi(argv[1]);
  char* mesh_name = argv[2];

  char geom_name[BUFSIZ];
  strncpy(geom_name,mesh_name,BUFSIZ);
  strncat(geom_name,".re2",BUFSIZ);

  char conn_name[BUFSIZ];
  strncpy(conn_name,mesh_name,BUFSIZ);
  strncat(conn_name,".co2",BUFSIZ);

  int color = MPI_UNDEFINED;
  if(rank<n_read) color=1;
  MPI_Comm comm_read;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm_read);

  long long *vl=NULL;
  double *coord=NULL;
  int nelt=0,ndim,nv;

  /* Read mesh data */
  if(color==1){
    struct comm comm; comm_init(&comm, comm_read);

    Mesh mesh;
    read_geometry(&mesh,geom_name,&comm);
    read_connectivity(mesh,conn_name,&comm);
    get_vertex_ids(&vl,mesh);
    get_vertex_coordinates(&coord,mesh);

    ndim = get_mesh_dim(mesh);
    nelt = get_mesh_nel(mesh);

    mesh_free(mesh);
    comm_free(&comm);
  }

  MPI_Bcast(&ndim,1,MPI_INT,0,MPI_COMM_WORLD);
  nv=(ndim==3)?8:4;

  /* Partition the mesh */
  int options[3];
  options[0] = 1; /* use custom options */
  options[1] = 2; /* debug level        */
  options[2] = 0; /* not used           */

  int *part=(int*)calloc(nelt,sizeof(int));
  int ierr=parRSB_partMesh(part,vl,coord,nelt,nv,options,MPI_COMM_WORLD);
  if(ierr){
    if(vl) free(vl);
    if(coord) free(coord);
    if(part) free(part);
    EXIT_ERROR();
  }

  /* Redistribute data */
  struct array elements;
  array_init(elm_data,&elements,nelt);
  elm_data *data;
  int e,n;
  for(data=elements.ptr,e=0; e<nelt; ++e) {
    data[e].proc=part[e];
    for(n=0; n<nv; ++n) {
      data[e].vtx[n] = vl[e*nv+n];
    }
  }
  elements.n=nelt;

  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);
  struct crystal cr; crystal_init(&cr, &comm);
  sarray_transfer(elm_data,&elements,proc,0,&cr);
  crystal_free(&cr);
  comm_free(&comm);

  nelt=elements.n;
  vl=(long long*) malloc(nv*nelt*sizeof(long long));
  for(data=elements.ptr, e=0; e<nelt; ++e) {
    for(n=0; n<nv; ++n){
      vl[e*nv+n]=data[e].vtx[n];
    }
  }
  array_free(&elements);

  printPartStat(vl, nelt, nv, MPI_COMM_WORLD);

  if(part) free(part);
  if(vl) free(vl);
  if(coord) free(coord);
}
