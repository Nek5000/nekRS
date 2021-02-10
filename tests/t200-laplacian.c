#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  struct comm comm;
  comm_init(&comm,MPI_COMM_WORLD);
  uint rank=comm.id,size=comm.np;

  if(argc!=2){
    if(rank==0) printf("Usage: ./%s <co2 file>\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  read_co2_mesh(&mesh,argv[1],&comm);

  genmap_handle gh; genmap_init(&gh,MPI_COMM_WORLD);
  GenmapSetNLocalElements(gh,mesh->nelt);
  GenmapSetNVertices(gh,mesh->nVertex);

  GenmapElements e=GenmapGetElements(gh);
  Element me      =MeshGetElements(mesh);
  GenmapInt i,j;
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i].vertex[j].globalId;

  GenmapVector weights,u,v;
  GenmapCreateVector(&weights,mesh->nelt);
  GenmapCreateVector(&u      ,mesh->nelt);
  GenmapCreateVector(&v      ,mesh->nelt);

  for(i=0;i<mesh->nelt;i++)
    u->data[i]=1.0;

  GenmapComm c=GenmapGetGlobalComm(gh);
  GenmapInitLaplacian(gh,c);
  GenmapLaplacian(gh,c,u,v);

  for(i=0;i<mesh->nelt;i++)
    assert(fabs(v->data[i])<GENMAP_TOL);

  GenmapDestroyVector(weights);
  GenmapDestroyVector(v);
  GenmapDestroyVector(u);

  mesh_free(mesh);
  genmap_finalize(gh);

  comm_free(&comm);
  MPI_Finalize();

  return 0;
}
