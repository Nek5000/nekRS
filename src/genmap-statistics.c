#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#include <genmap-impl.h>

#define MAXMETS 150
#define MAXLVLS  30
#define MAXSIZE (MAXMETS*MAXLVLS)

static double metrics[MAXMETS];
static double *stack;
static uint stack_size,stack_max;

void metric_init(){
  uint i; for(i=0; i<MAXMETS; i++)
    metrics[i]=0.0;
  stack=NULL;
  stack_size=stack_max=0;
}

void metric_finalize(){
  if(stack!=NULL)
    GenmapFree(stack);
}

void metric_acc(metric m,double count){ metrics[m]+=count; }

void metric_tic(struct comm *c,metric m){
  comm_barrier(c);
  metrics[m]-=comm_time();
}

void metric_toc(struct comm *c,metric m){
  metrics[m]+=comm_time();
  comm_barrier(c);
}

void metric_push_level(){
  assert(stack_size<=stack_max && "stack_size > stack_max");

  if(stack_size==stack_max){
    stack_max+=stack_size/2+1;
    GenmapRealloc(stack_max*MAXMETS,&stack);
  }

  uint i; for(i=0; i<MAXMETS; i++){
    stack[stack_size*MAXMETS+i]=metrics[i];
    metrics[i]=0.0;
  }
  stack_size++;
}

void metric_print(struct comm *c){
  double min[MAXSIZE],max[MAXSIZE],sum[MAXSIZE],buf[MAXSIZE];
  uint max_size=stack_size*MAXMETS;
  assert(max_size<=MAXSIZE);

  uint i; for(i=0; i<max_size; i++)
    min[i]=max[i]=sum[i]=stack[i];

  comm_allreduce(c,gs_double,gs_min,min,MAXSIZE,buf);// min
  comm_allreduce(c,gs_double,gs_max,max,MAXSIZE,buf);// max
  comm_allreduce(c,gs_double,gs_add,sum,MAXSIZE,buf);// sum

  for(i=0; i<max_size; i++)
    sum[i]/=c->np;

#define SUMMARY(i,m) sum[i*MAXMETS+m],min[i*MAXMETS+m],max[i*MAXMETS+m]

  for(i=0; i<stack_size; i++){
    if(c->id==0){
      printf("level=%02d\n",i);
      printf("  BINN1        : %g/%g/%g\n",SUMMARY(i,BINN1       ));
      printf("  BINN2        : %g/%g/%g\n",SUMMARY(i,BINN2       ));
      printf("  AXISLEN      : %g/%g/%g\n",SUMMARY(i,AXISLEN     ));
      printf("  LOCALSORT    : %g/%g/%g\n",SUMMARY(i,LOCALSORT   ));
      printf("  SETPROC      : %g/%g/%g\n",SUMMARY(i,SETPROC     ));
      printf("  RCBTRANSFER  : %g/%g/%g\n",SUMMARY(i,RCBTRANSFER ));
      printf("  COMMSPLIT    : %g/%g/%g\n",SUMMARY(i,COMMSPLIT));
      printf("  LOADBALANCE0 : %g/%g/%g\n",SUMMARY(i,LOADBALANCE0));
      printf("  LOADBALANCE1 : %g/%g/%g\n",SUMMARY(i,LOADBALANCE1));
      printf("  PARSORT      : %g/%g/%g\n",SUMMARY(i,PARSORT     ));
      printf("  UPDATEPROBE  : %g/%g/%g\n",SUMMARY(i,UPDATEPROBE ));
    }
  }
}

#undef SUMMARY

#undef MAXMETS
#undef MAXLVLS
#undef MAXSIZE
