#include <math.h>
#include <string.h>

#include <genmap-impl.h>
#include <gencon-impl.h>

#define readT(coords,buf,T,nVertex) do{\
  memcpy((coords),buf,sizeof(T)*nVertex);\
} while(0)

#define writeInt(dest,val) do{\
  memcpy(dest,&(val),sizeof(int));\
} while(0)

int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES]={0,1,3,2,4,5,7,6};
int PRE_TO_SYM_FACE[GC_MAX_FACES]={2,1,3,0,4,5};

int transferBoundaryFaces(Mesh mesh,struct comm *c){
  uint size=c->np;

  struct array *boundary=&mesh->boundary;
  BoundaryFace ptr=boundary->ptr;
  int nFaces=boundary->n;

  slong nelgt=mesh->nelgt;
  sint nelt=nelgt/size;
  sint nrem=nelgt-nelt*size;
  slong N=(size-nrem)*nelt;

  sint i; slong eid;
  for(i=0;i<nFaces;i++){
    eid=ptr[i].elementId;
    if(eid<N) ptr[i].proc=eid/nelt;
    else ptr[i].proc=(eid-N)/(nelt+1)+size-nrem;
  }

  struct crystal cr; crystal_init(&cr,c);
  sarray_transfer(struct Boundary_private,boundary,proc,1,&cr);
  crystal_free(&cr);
}

int readRe2Header(Mesh *mesh_,MPI_File file,struct comm *c){
  char version[BUFSIZ];
  int nelgt,nDim,nelgv,err;
  MPI_Status st;

  uint rank=c->id;
  uint size=c->np;
  MPI_Comm comm=c->c;

  char *buf=(char *)calloc(GC_RE2_HEADER_LEN+1,sizeof(char));
  if(rank==0){
    err=MPI_File_read(file,buf,GC_RE2_HEADER_LEN,MPI_BYTE,&st);
    if(err) return 1;
  }
  MPI_Bcast(buf,GC_RE2_HEADER_LEN,MPI_BYTE,0,comm);
  if(rank==0){
#if defined(GENMAP_DEBUG)
    printf("header: %s\n",buf);
#endif
  }
  sscanf(buf,"%5s %d %d %d",version,&nelgt,&nDim,&nelgv);

  int nVertex=(nDim==2)?4:8;
  int nelt=nelgt/size,nrem=nelgt-nelt*size;
  nelt+=(rank>=(size-nrem) ? 1: 0);

  MeshInit(mesh_,nelt,nDim); Mesh mesh=*mesh_;
  mesh->nelgt=nelgt;
  mesh->nelgv=nelgv;

  if(rank==0){
#if defined(GENMAP_DEBUG)
    printf("ndim,nvertex,nneighbors,nelgt,nelt:%d %d %d %lld %d\n",
      mesh->nDim,mesh->nVertex,mesh->nNeighbors,mesh->nelgt,mesh->nelt);
#endif
  }

  free(buf);

  return 0;
}

int readRe2Coordinates(Mesh mesh,MPI_File file,struct comm *c){
  uint rank=c->id;
  uint size=c->np;
  MPI_Comm comm=c->c;

  int nelt=mesh->nelt;
  int nelgt=mesh->nelgt;
  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;

  slong out[2][1],buff[2][1],in=nelt;
  comm_scan(out,c,gs_long,gs_add,&in,1,buff);
  slong start=out[0][0];

  int elemDataSize=nVertex*nDim*sizeof(double)+sizeof(double);
  int headerSize=GC_RE2_HEADER_LEN+sizeof(float);

  /* calculate read size for element data on each MPI rank */
  int readSize=nelt*elemDataSize;
  if(rank==0) readSize+=headerSize;

  char *buf=(char*)calloc(readSize,sizeof(char)),*buf0=buf;
  MPI_Status st;
  int err=MPI_File_read_ordered(file,buf,readSize,MPI_BYTE,&st);
  if(err) return 1;

  if(rank==0) buf0+=headerSize;

  /* initialize array */
  uint nUnits=nelt*nVertex;
  array_init(struct Point_private,&mesh->elements,nUnits);
  Point ptr=mesh->elements.ptr;

  /* read elements for each rank */
  double x[GC_MAX_VERTICES],y[GC_MAX_VERTICES],z[GC_MAX_VERTICES];
  int i,j,k;
  for(i=0;i<nelt;i++){
    // skip group id
    buf0+=sizeof(double);
    readT(x,buf0,double,nVertex); buf0+=sizeof(double)*nVertex;
    readT(y,buf0,double,nVertex); buf0+=sizeof(double)*nVertex;
    if(nDim==3){
      readT(z,buf0,double,nVertex);
      buf0+=sizeof(double)*nVertex;
    }

    for(k=0;k<nVertex;k++){
      j=PRE_TO_SYM_VERTEX[k];
      ptr->x[0]=x[j],ptr->x[1]=y[j];
      if(nDim==3)
        ptr->x[2]=z[j];
      ptr->elementId =start+i;
      ptr->sequenceId=nVertex*(start+i)+k;
      ptr->origin    =rank;
      ptr++;
    }
  }
  mesh->elements.n=nUnits;

#if defined(GENMAP_DEBUG)
  printf("io: rank=%d npts=%u\n",rank,nUnits);
#endif

  free(buf);
  return 0;
}

int readRe2Boundaries(Mesh mesh,MPI_File file,struct comm *c){
  uint rank=c->id;
  uint size=c->np;
  MPI_Comm comm=c->c;

  int nelt=mesh->nelt;
  int nelgt=mesh->nelgt;
  int nDim=mesh->nDim;
  int nVertex=mesh->nVertex;

  int elemDataSize=nVertex*nDim*sizeof(double)+sizeof(double);
  int headerSize=GC_RE2_HEADER_LEN+sizeof(float);

  MPI_Status st;
  char bufL[8];

  /* calculate offset for the curve side data */
  MPI_Offset curveOffset=headerSize+nelgt*elemDataSize;
  if(rank==0)
    MPI_File_read_at(file,curveOffset,bufL,sizeof(long),\
        MPI_BYTE,&st);
  MPI_Bcast(bufL,sizeof(long),MPI_BYTE,0,comm);

  double ncurvesD; readT(&ncurvesD,bufL,long,1);
  long ncurves=ncurvesD;

  /* calculate offset for boundary conditions data */
  MPI_Offset boundaryOffset=curveOffset+sizeof(long)+\
    sizeof(long)*8*ncurves;
  if(rank==0)
    MPI_File_read_at(file,boundaryOffset,bufL,sizeof(long),\
      MPI_BYTE,&st);
  MPI_Bcast(bufL,sizeof(long),MPI_BYTE,0,comm);

  double nbcsD; readT(&nbcsD,bufL,long,1); long nbcs=nbcsD;
  
  int nbcsLocal=nbcs/size,nrem=nbcs-nbcsLocal*size;
  nbcsLocal+=(rank>=(size-nrem) ? 1 : 0);

  slong out[2][1],buff[2][1],in=nbcsLocal;
  comm_scan(out,c,gs_long,gs_add,&in,1,buff);
  slong start=out[0][0];

  int offset=boundaryOffset+sizeof(long)+start*8*sizeof(long);
  int readSize=nbcsLocal*sizeof(long)*8;
  char *buf=calloc(readSize,sizeof(char)),*buf0=buf;
  MPI_File_read_at_all(file,offset,buf,readSize,MPI_BYTE,&st);

  double tmp[5];
  char cbc[4];
  struct Boundary_private boundary;
  sint i;
  for(i=0;i<nbcsLocal;i++){
    readT(tmp,buf0,long,1);buf0+=sizeof(long);
    boundary.elementId=tmp[0]-1;

    readT(tmp,buf0,long,1);buf0+=sizeof(long);
    boundary.faceId=PRE_TO_SYM_FACE[(long)tmp[0]-1];

    readT(tmp,buf0,long,5);buf0+=5*sizeof(long);
    readT(cbc,buf0,char,3);buf0+=sizeof(long);
    cbc[3]='\0';

    if(strcmp(cbc,GC_PERIODIC)==0){
      boundary.bc[0]=(long)tmp[0]-1;
      boundary.bc[1]=PRE_TO_SYM_FACE[(long)tmp[1]-1];
      array_cat(struct Boundary_private,&mesh->boundary,&boundary,1);
    }
  }
  free(buf);
}

int read_re2_mesh(Mesh *mesh_,char *fileName,struct comm *c){
  int nelt,nDim,nVertex;
  int errs=0;

  uint rank=c->id;
  uint size=c->np;
  MPI_Comm comm=c->c;

  MPI_File file;
  int err=MPI_File_open(comm,fileName,MPI_MODE_RDONLY,\
    MPI_INFO_NULL,&file);
  if(err){
    if(rank==0)
      fprintf(stderr,"%s:%d Error opening file: %s\n",
        __FILE__,__LINE__,fileName);
    errs++;
    MPI_Abort(comm,911);
  }

  readRe2Header(mesh_,file,c);
  Mesh mesh=*mesh_;
  readRe2Coordinates(mesh,file,c);
  readRe2Boundaries(mesh,file,c);
  transferBoundaryFaces(mesh,c);

  err=MPI_File_close(&file);
  if(err) errs++;

  MPI_Barrier(comm);

  return errs;
}

int read_co2_mesh(Mesh *mesh_,char *fname,struct comm *c){
  comm_ext comm=c->c;
  int size=c->np;
  int rank=c->id;

  MPI_File file;
  int err=MPI_File_open(comm,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,&file);
  if(err){
    if(rank==0)
      fprintf(stderr,"%s:%d Error opening file: %s for reading.\n",
        __FILE__,__LINE__,fname);
    MPI_Abort(comm,911);
  }

  char *buf=(char*)calloc(GC_CO2_HEADER_LEN+1,sizeof(char));
  MPI_Status st;
  if(rank==0){
    err=MPI_File_read(file,buf,GC_CO2_HEADER_LEN,MPI_BYTE,&st);
    if(err) return 1;
  }

  MPI_Bcast(buf,GC_CO2_HEADER_LEN,MPI_BYTE,0,comm);

  int nelgt,nelgv,nVertex;
  char version[6];
  sscanf(buf,"%5s%12d%12d%12d",version,&nelgt,&nelgv,&nVertex);

#if defined(GENMPA_DEBUG)
  if(rank==0)
    printf("%s %d %d %d\n",version,nelgt,nelgv,nVertex);
#endif

  //TODO: Assert version
  int nelt=nelgt/size,nrem=nelgt-nelt*size;
  nelt+=(rank>=(size-nrem) ? 1: 0);

  int nDim=(nVertex==8)?3:2;

  // Initialize the mesh structure
  MeshInit(mesh_,nelt,nDim); Mesh mesh=*mesh_;
  mesh->nelgt=nelgt;
  mesh->nelgv=nelgv;

#if defined(GENMPA_DEBUG)
  if(rank==0)
    printf("ndim/nvertex/nelgt/nelgv/nelt: %d/%d/%d/%d/%d\n",
      nDim,nVertex,nelgt,nelgv,nelt);
#endif

  slong out[2][1],buff[2][1],in[1];
  in[0]=nelt;
  comm_scan(out,c,gs_long,gs_add,&in,1,buff);
  slong start=out[0][0];

  int readSize=nelt*(nVertex+1)*sizeof(int);
  int headerSize=GC_CO2_HEADER_LEN+sizeof(float);
  if(rank==0)
    readSize+=headerSize;

  buf=(char*)realloc(buf,readSize*sizeof(char));
  err=MPI_File_read_ordered(file,buf,readSize,MPI_BYTE,&st);
  err=MPI_File_close(&file);

  char *buf0=buf;
  if(rank==0) buf0+=headerSize;

  Element ptr=mesh->elements.ptr;
  int i,j,tmp1,tmp2;
  for(i=0;i<nelt;i++){
    readT(&tmp1,buf0,int,1); buf0+=sizeof(int);
    for(j=0;j<nVertex;j++){
      ptr[i].vertex[j].elementId=tmp1;
      readT(&tmp2,buf0,int,1); buf0+=sizeof(int);
      ptr[i].vertex[j].globalId =tmp2;
    }
  }
  mesh->elements.n=nelt*nVertex;

  free(buf);

  return 0;
}

int read_co2_file(Mesh mesh,char *fname,struct comm *c){
  comm_ext comm=c->c;
  int rank=c->id;
  int size=c->np;

  MPI_File file;
  int err=MPI_File_open(comm,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,&file);
  if(err){
    if(rank==0)
      fprintf(stderr,"%s:%d Error opening file: %s for reading.\n",
        __FILE__,__LINE__,fname);
    MPI_Abort(comm,911);
  }

  char *buf=(char*)calloc(GC_CO2_HEADER_LEN+1,sizeof(char));
  MPI_Status st;

  if(rank==0){
    err=MPI_File_read(file,buf,GC_CO2_HEADER_LEN,MPI_BYTE,&st);
    if(err) return 1;
  }

  MPI_Bcast(buf,GC_CO2_HEADER_LEN,MPI_BYTE,0,comm);

  int nelgt,nelgv,nVertex;
  char version[6];
  sscanf(buf,"%5s%12d%12d%12d",version,&nelgt,&nelgv,&nVertex);

#if defined(GENMPA_DEBUG)
  if(rank==0)
    printf("%s %d %d %d\n",version,nelgt,nelgv,nVertex);
#endif

  //TODO: Assert version
  int nelt=nelgt/size,nrem=nelgt-nelt*size;
  nelt+=(rank<nrem ? 1: 0);

  int nDim=(nVertex==8)?3:2;

  // Initialize the mesh structure
  assert(mesh->nDim==nDim);
  assert(mesh->nVertex==nVertex);
  assert(mesh->nNeighbors==nDim);
  assert(mesh->nelgt==nelgt);
  assert(mesh->nelgv==nelgv);
  assert(mesh->nelt==nelt);
  assert(mesh->elements.n==nelt*nVertex);

#if defined(GENMPA_DEBUG)
  if(rank==0)
    printf("ndim/nvertex/nelgt/nelgv/nelt: %d/%d/%d/%d/%d\n",
      nDim,nVertex,nelgt,nelgv,nelt);
#endif

  slong out[2][1],buff[2][1],in[1];
  in[0]=nelt;
  comm_scan(out,c,gs_long,gs_add,&in,1,buff);
  slong start=out[0][0];

  int readSize=nelt*(nVertex+1)*sizeof(int);
  int headerSize=GC_CO2_HEADER_LEN+sizeof(float);
  if(rank==0)
    readSize+=headerSize;

  buf=(char*)realloc(buf,readSize*sizeof(char));
  err=MPI_File_read_ordered(file,buf,readSize,MPI_BYTE,&st);
  err=MPI_File_close(&file);

  char *buf0=buf;
  if(rank==0) buf0+=headerSize;

  Point ptr=mesh->elements.ptr;
  int i,j,tmp1,tmp2;
  for(i=0;i<nelt;i++){
    readT(&tmp1,buf0,int,1); buf0+=sizeof(int);
    for(j=0;j<nVertex;j++){
      ptr[i*nVertex+j].elementId=tmp1;
      readT(&tmp2,buf0,int,1); buf0+=sizeof(int);
      ptr[i*nVertex+j].globalId =tmp2;
    }
  }

  free(buf);

  return 0;
}

int write_co2_file(Mesh mesh,char *fileName,struct comm *c){
  const char version[5]="#v001";
  const float test=6.54321;

  uint rank=c->id;
  uint size=c->np;
  MPI_Comm comm=c->c;

  int nVertex=mesh->nVertex;
  int nDim=mesh->nDim;
  sint nelt=mesh->nelt;
  slong nelgt=mesh->nelgt;
  slong nelgv=mesh->nelgv;

  int errs=0;

  MPI_File file;
  int err=MPI_File_open(comm,fileName,\
    MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&file);
  if(err){
    if(rank==0)
      fprintf(stderr,"%s:%d Error opening file: %s for writing.\n",
        __FILE__,__LINE__,fileName);
    errs++;
    MPI_Abort(comm,911);
  }

  slong out[2][1],buff[2][1],in=nelt;
  comm_scan(out,c,gs_long,gs_add,&in,1,buff);
  slong start=out[0][0];

  int writeSize=nelt*(nVertex+1)*sizeof(int);
  int headerSize=GC_CO2_HEADER_LEN+sizeof(float);
  if(rank==0) writeSize+=headerSize;

  char *buf=(char*)calloc(writeSize,sizeof(char)),*buf0=buf;
  MPI_Status st;
  if(rank==0){
    sprintf(buf0,"%5s%12d%12d%12d",version,(int)nelgt,\
        (int)nelgv,nVertex);
#if defined(GENMAP_DEBUG)
    printf("%5s%12d%12d%12d\n",version,(int)nelgt,(int)nelgv,nVertex);
#endif
    memset(buf0+strlen(buf0),' ',GC_CO2_HEADER_LEN-strlen(buf0));
    buf0[GC_CO2_HEADER_LEN]='\0';
    buf0+=GC_CO2_HEADER_LEN;
    memcpy(buf0,&test,sizeof(float)),buf0+=sizeof(float);
  }

  Point ptr=mesh->elements.ptr;
  int i,k,temp;
  for(i=0;i<nelt;i++){
    temp=ptr->elementId+1;
    writeInt(buf0,temp); buf0+=sizeof(int);
#if defined(GENMAP_DEBUG)
    printf("%d ",temp);
#endif
    for(k=0;k<nVertex;k++){
      temp=ptr->globalId+1;
      writeInt(buf0,temp); buf0+=sizeof(int);
#if defined(GENMAP_DEBUG)
      printf("%d ",temp);
#endif
      ptr++;
    }
#if defined(GENMAP_DEBUG)
    printf("\n");
#endif
  }

  err=MPI_File_write_ordered(file,buf,writeSize,MPI_BYTE,&st);
  if(err) errs++;

  err=MPI_File_close(&file);
  if(err) errs++;

  MPI_Barrier(comm);
  free(buf);

  return errs;
}

#undef readT
#undef writeInt
