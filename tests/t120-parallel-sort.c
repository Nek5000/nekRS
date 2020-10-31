#include <genmap-impl.h>
#include <sort-impl.h>
#include <time.h>

typedef struct{
  double ds;
  slong dl;
  uint proc;
} Data;

#define N 10

void check(struct array *arr,struct comm *c){
  uint size=c->np,rank=c->id;

  Data *ptr=arr->ptr;
  uint n=arr->n,i;
  for(i=0; i<n-1; i++){
    assert(ptr[i].ds<=ptr[i+1].ds && "Field ds is not sorted");
    if(i%2==0)
      assert(ptr[i].dl<=ptr[i+1].dl && "Field dl is not sorted");
  }

  uint root=0;

  struct array a; array_init(Data,&a,2); Data *p=a.ptr;
  p[0]     =ptr[0],p[1]     =ptr[n-1];
  p[0].proc=root  ,p[1].proc=root;
  p[0].dl  =2*rank,p[1].dl  =2*rank+1;
  a.n=2;

  struct crystal cr; crystal_init(&cr,c);

  sarray_transfer(Data,&a,proc,0,&cr);
  if(rank==root)
    assert(a.n==2*size);

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(Data,a.ptr,a.n,dl,1,&buf);
  buffer_free(&buf);

  ptr=a.ptr;
  if(rank==root)
    for(i=0; i<a.n-1; i++)
      assert(ptr[i].ds<=ptr[i+1].ds && "Field ds is not sorted globally");

  crystal_free(&cr);
  array_free(&a);
}

int main(int argc,char *argv[]){
  MPI_Init(&argc,&argv);
  struct comm c; comm_init(&c,MPI_COMM_WORLD);

  struct array arr; array_init(Data,&arr,N);

  srand(time(0));

  int i,cnt; Data d; Data *ptr=arr.ptr;
  for(i=cnt=0; i<N/2; i++){
      d.dl=rand()%100,d.ds=(rand()%100)/100.0,ptr[cnt++]=d;
      d.dl=rand()%100,ptr[cnt++]=d;
  }
  arr.n=N;

  parallel_sort_2(Data,&arr,ds,gs_double,dl,gs_long,0,1,&c);
  check(&arr,&c);

  //TODO: Test hypercube sort

  array_free(&arr);
  comm_free(&c);

  MPI_Finalize();

  return 0;
}
