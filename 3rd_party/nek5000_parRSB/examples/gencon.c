#include <stdio.h>
#include <stdlib.h>

#include <gencon-impl.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);

  if(argc<2){
    if(comm.id==0)
      printf("Usage: ./%s foo.re2 [tol]\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  read_geometry(&mesh,argv[1],&comm);

  findMinNeighborDistance(mesh);

  double tol=(argc>2)?atof(argv[2]):0.2;

  findSegments(mesh,&comm,tol,0);
  setGlobalID(mesh,&comm);
  sendBack(mesh,&comm);
  matchPeriodicFaces(mesh,&comm);

  char co2FileName[BUFSIZ]; strncpy(co2FileName,argv[1],BUFSIZ);
  int len=strlen(co2FileName); assert(len>4 && len<BUFSIZ);
  co2FileName[len-2]='o',co2FileName[len-3]='c';

  write_connectivity(mesh,co2FileName,&comm);

  mesh_free(mesh);
  comm_free(&comm);
  MPI_Finalize();

  return 0;
}
