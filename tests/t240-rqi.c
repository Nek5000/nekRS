#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <genmap-impl.h>
#include <gencon-impl.h>
#include <genmap-multigrid-precon.h>
#include <parRSB.h>

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);

  struct comm comm; comm_init(&comm,MPI_COMM_WORLD);
  int rank=comm.id,size=comm.np;

  buffer buf; buffer_init(&buf,1024);

  if(argc<3){
    if(rank==0)
      printf("Usage: ./%s <re2 file> <co2 file> global local\n",argv[0]);
    MPI_Finalize();
    exit(1);
  }

  Mesh mesh;
  read_re2_mesh(&mesh,argv[1],&comm);
  read_co2_file( mesh,argv[2],&comm);

  GenmapInt i,j;

  int rcb_g=(argc>3)?atoi(argv[3]):1;
  int rcb_l=(argc>4)?atoi(argv[4]):1;

  //partition
  int *part; GenmapMalloc(mesh->nelt,&part);
  int *seq; GenmapMalloc(mesh->nelt,&seq);
  uint *upart; GenmapMalloc(mesh->nelt*mesh->nVertex,&upart);
  double *coords;
  GenmapMalloc(mesh->nelt*mesh->nVertex*mesh->nDim,&coords);

  int nDim=mesh->nDim;
  Point me=(Point)MeshGetElements(mesh);

  if(rcb_g){
    int cnt=0;
    for(i=0; i<mesh->nelt; i++){
      for(j=0; j<mesh->nVertex; j++){
        coords[cnt++]=me[i*mesh->nVertex+j].x[0];
        coords[cnt++]=me[i*mesh->nVertex+j].x[1];
        if(nDim==3)
          coords[cnt++]=me[i*mesh->nVertex+j].x[2];
      }
    }

    int options[3]; options[0]=options[1]=options[2]=0;
    parRCB_partMesh(part,seq,coords,mesh->nelt,mesh->nVertex,options,
      MPI_COMM_WORLD);

    for(i=0; i<mesh->nelt; i++)
      for(j=0; j<mesh->nVertex; j++)
        upart[i*mesh->nVertex+j]=part[i];

    struct crystal cr; crystal_init(&cr,&comm);
    sarray_transfer_ext(struct Point_private,&mesh->elements,upart,
        sizeof(uint),&cr);
    crystal_free(&cr);

    sarray_sort(struct Point_private,mesh->elements.ptr,
        mesh->elements.n,sequenceId,1,&buf);
    mesh->nelt=mesh->elements.n/mesh->nVertex;
  }

  if(rcb_l){
    struct array a; array_init(elm_rcb,&a,mesh->nelt); a.n=mesh->nelt;
    elm_rcb *ptr=a.ptr;
    for(i=0; i<mesh->nelt; i++){
      ptr[i].orig=i;
      ptr[i].coord[0]=0.0;
      ptr[i].coord[1]=0.0;
      if(nDim==3)
        ptr[i].coord[2]=0.0;
      for(j=0; j<mesh->nVertex; j++){
        ptr[i].coord[0]+=me[i*mesh->nVertex+j].x[0];
        ptr[i].coord[1]+=me[i*mesh->nVertex+j].x[1];
        if(nDim==3)
          ptr[i].coord[2]+=me[i*mesh->nVertex+j].x[2];
      }
      ptr[i].coord[0]/=mesh->nVertex;
      ptr[i].coord[1]/=mesh->nVertex;
      if(nDim==3)
        ptr[i].coord[2]/=mesh->nVertex;
    }

    uint s1=0,e1=mesh->nelt;
    rcb_local(&a,s1,e1,mesh->nDim,&buf);
    ptr=a.ptr;

    for(i=0; i<mesh->nelt; i++) ptr[i].proc=i;

    sarray_sort(elm_rcb,a.ptr,a.n,orig,0,&buf);
    ptr=a.ptr;

    Point pp=mesh->elements.ptr;
    int cnt=0;
    for(i=0; i<mesh->nelt; i++)
      for(j=0; j<mesh->nVertex; j++){
        pp[cnt].proc=ptr[i].proc;
        cnt++;
      }
    array_free(&a);

    sarray_sort_2(struct Point_private,mesh->elements.ptr,
        mesh->elements.n,proc,0,sequenceId,1,&buf);
  }

  free(upart);
  free(part);
  free(seq);
  free(coords);
  buffer_free(&buf);

  GenmapHandle gh; GenmapInit(&gh,MPI_COMM_WORLD);
  GenmapSetNLocalElements(gh,mesh->nelt);
  GenmapSetNVertices(gh,mesh->nVertex);

  /* Setup mesh */
  GenmapElements e=GenmapGetElements(gh);
  me=(Point)MeshGetElements(mesh);
  for(i=0;i<mesh->nelt;i++)
    for(j=0;j<mesh->nVertex;j++)
      e[i].vertices[j]=me[i*mesh->nVertex+j].globalId;

  GenmapComm c=GenmapGetGlobalComm(gh);
  GenmapInitLaplacian(gh,c);

  GenmapVector x; GenmapCreateVector(&x,mesh->nelt);
  GenmapVector r; GenmapCreateVector(&r,mesh->nelt);

  srand(time(0));
  for(i=0; i<mesh->nelt; i++){
    x->data[i]=rand()%100/50.;
  }

  GenmapLong nelg=GenmapGetNGlobalElements(gh);
  GenmapOrthogonalizebyOneVector(gh,c,x,nelg);

  GenmapScalar norm=GenmapDotVector(x,x);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);

  GenmapScalar normi=1.0/sqrt(norm);
  GenmapAxpbyVector(x,x,0.0,x,normi);

  norm=GenmapDotVector(x,x);
  GenmapGop(c,&norm,1,GENMAP_SCALAR,GENMAP_SUM);

  mgData d; mgSetup(c,c->M,&d); d->h=gh;

  rqi(gh,c,d,x,30,1,r);

  mgFree(d);

  GenmapDestroyVector(x); GenmapDestroyVector(r);

  GenmapFinalize(gh);
  MeshFree(mesh);

  comm_free(&comm);
  MPI_Finalize();

  return 0;
}
