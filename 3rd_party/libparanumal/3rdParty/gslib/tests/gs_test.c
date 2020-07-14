#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"
#include "mem.h"
#include "gs_defs.h"
#include "gs.h"

typedef double T;
const gs_dom dom = gs_double;

static void test(const struct comm *comm, gs_method method)
{
  struct gs_data *gsh;
  const uint np = comm->np;
  slong *id = tmalloc(slong,np+4);
  T *v = tmalloc(T,np+4);
  uint i;
  id[0] = -(slong)(np+10+3*comm->id);
  for(i=0;i<np;++i) id[i+1] = -(sint)(i+1);
  id[np+1] = comm->id+1;
  id[np+2] = comm->id+1;
  id[np+3] = np-comm->id;
  gsh = gs_setup(id,np+4,comm,0,method,1);
  free(id);

  /* non-blocking api - original test */
  if(comm->id==0) printf("\nTesting non-blocking api ...\n");
  for(i=0;i<np+4;++i) v[i] = 1;
  int handle;
  igs(v,dom,gs_add,0,gsh,0,&handle);
  gs_wait (handle);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);
  if(comm->id==0) printf("\n");

  for(i=0;i<np+4;++i) v[i] = 1;
  igs(v,dom,gs_add,1,gsh,0,&handle);
  gs_wait (handle);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);

  /* blocking api - original test */
  if(comm->id==0) printf("\nTesting blocking api ...\n");
  for(i=0;i<np+4;++i) v[i] = 1;
  gs(v,dom,gs_add,0,gsh,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);
  if(comm->id==0) printf("\n");

  for(i=0;i<np+4;++i) v[i] = 1;
  gs(v,dom,gs_add,1,gsh,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);
  gs_free(gsh);
  free(v);

  /* non-blocking gs, gs_many and gs_vec */
  uint neighbors = 3; //max neighbors, counting ownself
  slong id1[neighbors];
  slong me = (slong)comm->id;
  uint count=0;
  if(me>0) id1[count++]=me;
  id1[count++]=me+1;
  if(me<np-1) id1[count++]=me+2;
  neighbors=count;

  struct gs_data *gsh1;
  gsh1 = gs_setup(id1,neighbors,comm,0,method,0);

  T answer[np];
  if   (np==1) answer[0]=1.0;
  else {
    answer[0]=answer[np-1]=2.0;
    for(i=1;i<np-1;i++) answer[i]=3.0;
  }

  T u[neighbors];
  for(i=0;i<neighbors;i++) u[i]=1.0;

  igs(u,dom,gs_add,0,gsh,0,&handle);
  gs_wait(handle);

  for(i=0;i<neighbors;i++) assert(fabs(u[i]-answer[id1[i]-1])<1e-16);

  T w[neighbors][2];
  for(i=0;i<neighbors;i++) w[i][0]=1.0,w[i][1]=2.0;
  igs_vec(w,2,dom,gs_add,0,gsh,0,&handle);
  gs_wait(handle);

  for(i=0;i<neighbors;i++) assert(fabs(w[i][0]-answer[id1[i]-1])<1e-16);
  for(i=0;i<neighbors;i++) assert(fabs(w[i][1]-2*answer[id1[i]-1])<1e-16);

  T x1[neighbors], x2[neighbors];
  for(i=0;i<neighbors;i++) x1[i]=1.0,x2[i]=2.0;
  T *x[2] = {x1, x2};
  igs_many((void*)x,2,dom,gs_add,0,gsh,0,&handle);
  gs_wait(handle);

  for(i=0;i<neighbors;i++) assert(fabs(x1[i]-answer[id1[i]-1])<1e-16);
  for(i=0;i<neighbors;i++) assert(fabs(x2[i]-2*answer[id1[i]-1])<1e-16);
}

int main(int narg, char *arg[])
{
  comm_ext world; int np;
  struct comm comm;

#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif

  comm_init(&comm,world);

  test(&comm,gs_all_reduce);
  test(&comm,gs_pairwise);

  comm_free(&comm);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
