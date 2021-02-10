#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "gslib.h"
#include "rand_elt_test.h"

#define D 3

#if D==3
#define INITD(a,b,c) {a,b,c}
#define MULD(a,b,c) ((a)*(b)*(c))
#define INDEXD(a,na, b,nb, c) (((c)*(nb)+(b))*(na)+(a))
#define findpts_data    findpts_data_3
#define findpts_setup   findpts_setup_3
#define findpts_free    findpts_free_3
#define findpts         findpts_3
#define findpts_eval    findpts_eval_3
#define findptsms_setup findptsms_setup_3
#define findptsms_free  findptsms_free_3
#define findptsms       findptsms_3
#define findptsms_eval  findptsms_eval_3
#elif D==2
#define INITD(a,b,c) {a,b}
#define MULD(a,b,c) ((a)*(b))
#define INDEXD(a,na, b,nb, c) ((b)*(na)+(a))
#define findpts_data    findpts_data_2
#define findpts_setup   findpts_setup_2
#define findpts_free    findpts_free_2
#define findpts         findpts_2
#define findpts_eval    findpts_eval_2
#define findptsms_setup findptsms_setup_2
#define findptsms_free  findptsms_free_2
#define findptsms       findptsms_2
#define findptsms_eval  findptsms_eval_2
#endif


#define NR 5
#define NS 7
#define NT 6
#define K 4
#define NEL MULD(K,K,K)
#define TN 4

#define NPT_MAX 256
#define BBOX_TOL 0.01
#define NEWT_TOL 1024*DBL_EPSILON
#define LOC_HASH_SIZE NEL*MULD(NR,NS,NT)
#define GBL_HASH_SIZE NEL*MULD(NR,NS,NT)

/*
#define NPT_MAX 256
#define BBOX_TOL 1.00
#define NEWT_TOL 1024*DBL_EPSILON
#define LOC_HASH_SIZE NEL*3
#define GBL_HASH_SIZE NEL*3
*/

static uint np, id, idsess,np1,id1;
static double xdom0,xdom1;
static double volv,totalvolv;

static const unsigned nr[D] = INITD(NR,NS,NT);
static const unsigned mr[D] = INITD(2*NR,2*NS,2*NT);
static double zr[NR], zs[NS], zt[NT];
static double x3[D][MULD(3,3,3)];
static double mesh[D][NEL*MULD(NR,NS,NT)];
static double distfint[NEL*MULD(NR,NS,NT)];
static double gllsid[NEL*MULD(NR,NS,NT)];
static const double *const elx[D] = INITD(mesh[0],mesh[1],mesh[2]);
static const double *const dfint = {distfint};
static const double *const gllsi = {gllsid};

struct pt_data { double x[D], r[D], dist2, ex[D], esid, esidm; uint code, proc, el, sid, ptmarker, sidm;};
static struct array testp;

static struct crystal cr;
static struct crystal cr1;

static double quad_eval(const double coef[MULD(3,3,3)], const double r[D])
{
  double lr0[D], lr1[D], lr2[D];
  unsigned d;
  for(d=0;d<D;++d) lr0[d]=r[d]*(r[d]-1)/2,
                   lr1[d]=(1+r[d])*(1-r[d]),
                   lr2[d]=r[d]*(r[d]+1)/2;
  #define EVALR(base) ( coef [base   ]*lr0[0] \
                       +coef [base+ 1]*lr1[0] \
                       +coef [base+ 2]*lr2[0] )
  #define EVALS(base) ( EVALR(base   )*lr0[1] \
                       +EVALR(base+ 3)*lr1[1] \
                       +EVALR(base+ 6)*lr2[1] )
  #define EVALT(base) ( EVALS(base   )*lr0[2] \
                       +EVALS(base+ 9)*lr1[2] \
                       +EVALS(base+18)*lr2[2] )
  #if D==2
  #  define EVAL() EVALS(0)
  #elif D==3
  #  define EVAL() EVALT(0)
  #endif
  
  return EVAL();
  
  #undef EVAL
  #undef EVALT
  #undef EVALS
  #undef EVALR
}

static void test_mesh(void)
{
  const uint pn = ceil(pow(np1,1.0/D));
  const uint pi=id1%pn, pj=(id1/pn)%pn;
  #if D==3
  const uint pk=(id1/pn)/pn;
  #endif
  const double pfac = 1.0/pn;
  const double pbase[D] = INITD(-1+2*pfac*pi, -1+2*pfac*pj, -1+2*pfac*pk);
  const double fac = 1.0/K, step = 2.0/(K*(TN-1));
  unsigned ki,kj;
  #if D==3
  unsigned kk;
  #endif
  struct pt_data *out = testp.ptr;
  testp.n = NEL*MULD(TN,TN,TN);
  memset(testp.ptr,0,testp.n*sizeof(struct pt_data));
  uint idx;
  idx = 0;
  #if D==3
  for(kk=0;kk<K;++kk)
  #endif
  for(kj=0;kj<K;++kj) for(ki=0;ki<K;++ki) {
    unsigned i,j;
    double r[D], base[D] = INITD(-1+2*fac*ki,-1+2*fac*kj,-1+2*fac*kk);
    #if D==3
    unsigned k;
    for(k=0;k<TN;++k) { r[2] = pbase[2]+pfac*(1+base[2]+step*k);
    #endif
    for(j=0;j<TN;++j) { r[1] = pbase[1]+pfac*(1+base[1]+step*j);
    for(i=0;i<TN;++i) { r[0] = pbase[0]+pfac*(1+base[0]+step*i);
      out->proc = rand()%np1;
      out->x[0] = quad_eval(x3[0],r);
      out->x[1] = quad_eval(x3[1],r);
      #if D==3
      out->x[2] = quad_eval(x3[2],r);
      #endif

      out->sidm=0;
      if (idsess==0) {out->sidm=1;}
      out->sid = 2;
      out->ptmarker=2;

      if (out->x[0]<xdom0 && out->x[0]>xdom1)
      {
       if (idx%2==0) {out->ptmarker=1;out->sid=idsess;}
       else {out->ptmarker=0;out->sid=2;};
      }
      if (idsess==0 && fabs(out->x[0]-xdom0)<1.e-10) {out->ptmarker=1;out->sid=idsess;};
      if (idsess==1 && fabs(out->x[0]-xdom1)<1.e-10) {out->ptmarker=1;out->sid=idsess;};
      ++idx;
      ++out;
    }}
    #if D==3
    }
    #endif
  }
  sarray_transfer(struct pt_data,&testp,proc,1,&cr1);
  if(0)
    printf("%u: %u shuffled points\n",id,(unsigned)testp.n);
}

static void print_ptdata(const struct comm *const comm)
{
  uint notfound=0;
  uint notsid  =0;
  double dist2max=0, ed2max=0;
  const struct pt_data *pt = testp.ptr, *const end = pt+testp.n;
  for(;pt!=end;++pt) {
    if(0&&id==0)
    printf("code=%u, proc=%u, el=%u, dist2=%g, r=(%.17g,%.17g"
           #if D==3
           ",%.17g"
           #endif
           "), "
           "x=(%.17g,%.17g"
           #if D==3
           ",%.17g"
           #endif
           "), ex=(%.17g,%.17g"
           #if D==3
           ",%.17g"
           #endif
           ")\n",
      pt->code,pt->proc,pt->el,pt->dist2,
      pt->r[0],pt->r[1],
      #if D==3
      pt->r[2],
      #endif
      pt->x[0],pt->x[1],
      #if D==3
      pt->x[2],
      #endif
      pt->ex[0],pt->ex[1]
      #if D==3
      ,pt->ex[2]
      #endif
      );
    if(pt->code==2) ++notfound;
    else {
      double ed2=0, dx;
      unsigned d; for(d=0;d<D;++d) dx=pt->x[d]-pt->ex[d], ed2+=dx*dx;
      dist2max=pt->dist2>dist2max?pt->dist2:dist2max;
      ed2max=ed2>ed2max?ed2:ed2max;
    }
    uint corsid;
    corsid=0;
    if (pt->ptmarker==0) {
     corsid = 0;
     if (pt->x[0] > 0.5*(xdom0+xdom1)) corsid = 1;
     if (fabs(pt->x[0]-0.5*(xdom0+xdom1))>1.e-10 && (fabs(pt->esid-corsid)>1.e-2)) ++notsid;
    }
    else if (pt->ptmarker==1)
    {
    corsid=0;
    if (idsess==0) {corsid=1;}
    if ((pt->x[0]<xdom0 && pt->x[0]>xdom1) || fabs(pt->x[0]-xdom0)<1.e-10 || fabs(pt->x[0]-xdom1)<1.e-10) {
       if (fabs(pt->esid-corsid)>1.e-2) ++notsid;
        }
    else 
      {if (fabs(pt->esid-idsess)>1.e-2) ++notsid;
      }
    }
  }
  {
    double distmax=sqrt(dist2max), edmax=sqrt(ed2max);
    slong total=testp.n;
    if(0)
    printf("%u: maximum distance = %g (adv), %g (eval);"
           " %u/%u points not found in session %u\n",
         (unsigned)id, distmax, edmax,
         (unsigned)notfound, (unsigned)testp.n,idsess);
    distmax = comm_reduce_double(comm,gs_max,&distmax,1);
    edmax   = comm_reduce_double(comm,gs_max,&edmax  ,1);
    notfound = comm_reduce_sint(comm,gs_add,(sint*)&notfound,1);
    total    = comm_reduce_slong(comm,gs_add,&total,1);
    if(id1==0)
      printf("maximum distance = %g (adv), %g (eval);"
           " %u/%lu points not found in session %u\n",
           distmax, edmax,
           (unsigned)notfound, (unsigned long)total, idsess);
    if(id1==0)
      printf("Number of points identified with wrong session id: %u/%lu \n",
           (unsigned)notsid, (unsigned long)total);
  }
}

static void print_ptdata_match(const struct comm *const comm)
{
  uint notfound=0;
  const struct pt_data *pt = testp.ptr, *const end = pt+testp.n;
  for(;pt!=end;++pt) {
    if (idsess==0) {
      double diff = fabs(pt->esidm-1);
      if (pt->x[0]>xdom1 && (pt->code==2 || diff>1.e-2)) ++notfound;
    }
    else {
      double diff = fabs(pt->esidm-0);
      if (pt->x[0]<xdom0 && (pt->code==2 || diff>1.e-2)) ++notfound;
    }
    }
    slong total=testp.n;
    notfound = comm_reduce_sint(comm,gs_add,(sint*)&notfound,1);
    total    = comm_reduce_slong(comm,gs_add,&total,1);
    if(id1==0)
      printf("Number of points identified with wrong session id using findpts session matching: %u/%lu \n",
           (unsigned)notfound, (unsigned long)total);
}

static void rand_mesh1(double x0, double x1, double y0, double y1, double z0, double z1)
{
  const uint pn = ceil(pow(np1,1.0/D));
  const uint pi=id1%pn, pj=(id1/pn)%pn;
  #if D==3
  const uint pk=(id1/pn)/pn;
  #endif
  const double pfac = 1.0/pn;
  const double pfac0 = pfac;
  const double pfac1 = 1.0;
  const double pfac2 = 1.0;
  const double pbase[D] = INITD(-1+2*pfac*pi, -1+2*pfac*pj, -1+2*pfac*pk);
  const double fac = 1.0/K;
  const double z3[3] = {-1,0,1};
  unsigned ki,kj;
  unsigned kk;
  if(id==0) printf("Global division: %u^%d\n",(unsigned)pn,D);
  if(id==0) printf("Global division: %u^%d\n",(unsigned)pn,D);

  double dxl,dyl,dzl;
  dxl = x1-x0;
  dyl = y1-y0;
  dzl = z1-z0;

  unsigned ii,jj;
  int idx;
  idx = 0;
  #if D==3 
  for (ii=0;ii<3;++ii)
  #endif
  for (jj=0;jj<3;++jj)
  for (kk=0;kk<3;++kk)
  {
     x3[0][idx] = x0 + dxl*(z3[kk]+1.)/2.;
     x3[1][idx] = y0 + dyl*(z3[jj]+1.)/2.;
  #if D==3 
     x3[2][idx] = z0 + dzl*(z3[ii]+1.)/2.;
  #endif
   ++idx;
  }

  #if D==3
  for(kk=0;kk<K;++kk)
  #endif
  for(kj=0;kj<K;++kj) for(ki=0;ki<K;++ki) {
    unsigned off = INDEXD(ki,K, kj,K, kk)*MULD(NR,NS,NT);
    unsigned i,j;
    double r[D], base[D] = INITD(-1+2*fac*ki,-1+2*fac*kj,-1+2*fac*kk);
    #if D==3
    unsigned k;
    for(k=0;k<NT;++k) { r[2]=pbase[2]+pfac2*(1+base[2]+fac*(1+zt[k]));
    #endif
    for(j=0;j<NS;++j) { r[1]=pbase[1]+pfac1*(1+base[1]+fac*(1+zs[j]));
    for(i=0;i<NR;++i) { r[0]=pbase[0]+pfac0*(1+base[0]+fac*(1+zr[i]));
      idx = off+INDEXD(i,NR, j,NS, k);
      mesh[0][idx] = quad_eval(x3[0],r);
      mesh[1][idx] = quad_eval(x3[1],r);
      distfint[idx] = idsess==0?xdom0-mesh[0][idx]:mesh[0][idx]-xdom1;
      gllsid[idx] = idsess;
      #if D==3
      mesh[2][idx] = quad_eval(x3[2],r);
    }
    #endif
    }}
  }
}

static void test(const struct comm *const comm, const struct comm *const comm1)
{
  const double *x_base[D];
  const unsigned x_stride[D] = INITD(sizeof(struct pt_data),
                                     sizeof(struct pt_data),
                                     sizeof(struct pt_data));
  struct findpts_data *fd;
  struct pt_data *pt;
  unsigned d;
  if(id==0) printf("Initializing mesh\n");
  xdom0 = 1.0;
  xdom1 = 0.5;
  double xmn,xmx,ymn= 2.,ymx= 3.,zmn= 0.,zmx= 5.;
  if (idsess==0) {
    xmn=-3;xmx=xdom0;
  }
  else {
    xmn=xdom1;xmx=4.;
  }
  rand_mesh1(xmn,xmx,ymn,ymx,zmn,zmx);
  test_mesh();
  pt = testp.ptr;
  if(id==0) printf("calling findpts_setup\n");
  fd=findptsms_setup(comm,elx,nr,NEL,mr,BBOX_TOL,
                   LOC_HASH_SIZE,GBL_HASH_SIZE,
                   NPT_MAX,NEWT_TOL,&idsess,dfint);
  if(id==0) printf("calling findpts\n");
  x_base[0]=pt->x, x_base[1]=pt->x+1;
  #if D==3
  x_base[2]=pt->x+2;
  #endif
   uint match = 0;
   findptsms(&pt->code ,sizeof(struct pt_data),
             &pt->proc ,sizeof(struct pt_data),
             &pt->el   ,sizeof(struct pt_data),
             pt->r     ,sizeof(struct pt_data),
             &pt->dist2,sizeof(struct pt_data),
             x_base    ,x_stride, 
             &pt->sid,sizeof(struct pt_data),
             &match, testp.n, fd);
  for(d=0;d<D;++d) {
    if(id==0) printf("calling findpts_eval (%u)\n",d);
    findptsms_eval(&pt->ex[d],sizeof(struct pt_data),
                   &pt->code ,sizeof(struct pt_data),
                   &pt->proc ,sizeof(struct pt_data),
                   &pt->el   ,sizeof(struct pt_data),
                   pt->r     ,sizeof(struct pt_data),
                   testp.n   , mesh[d], fd);
  }
    findptsms_eval(&pt->esid   ,sizeof(struct pt_data),
                   &pt->code   ,sizeof(struct pt_data),
                   &pt->proc   ,sizeof(struct pt_data),
                   &pt->el     ,sizeof(struct pt_data),
                   pt->r       ,sizeof(struct pt_data),
                   testp.n     , gllsi, fd);
   print_ptdata(comm1);
   pt = testp.ptr;
// Test finding points in specified session ID
   match = 1;
   findptsms(&pt->code ,sizeof(struct pt_data),
             &pt->proc ,sizeof(struct pt_data),
             &pt->el   ,sizeof(struct pt_data),
             pt->r     ,sizeof(struct pt_data),
             &pt->dist2,sizeof(struct pt_data),
             x_base    ,x_stride,
             &pt->sidm ,sizeof(struct pt_data),
             &match, testp.n, fd);
   findptsms_eval(&pt->esidm  ,sizeof(struct pt_data),
                  &pt->code   ,sizeof(struct pt_data),
                  &pt->proc   ,sizeof(struct pt_data),
                  &pt->el     ,sizeof(struct pt_data),
                  pt->r       ,sizeof(struct pt_data),
                  testp.n, gllsi, fd);
  print_ptdata_match(comm1);
  findpts_free(fd);
}


int main(int narg, char *arg[])
{
  comm_ext world;
  struct comm comm;
  comm_ext world1;
  struct comm comm1;
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
#else
  world=0;
#endif

  comm_init(&comm,world);
  id=comm.id, np=comm.np;

  if (np<2) {
  printf("Please use atleast 2 processors\n");
  #ifdef MPI
    MPI_Finalize();
  #endif
  return 0;
  }
  else
  {
  idsess = 1;
  if (id%2==0) {idsess=0;}
  } 
  MPI_Comm_split(world,idsess,id,&world1);
  comm_init(&comm1,world1);
  id1 = comm1.id, np1=comm1.np;

  lobatto_nodes(zr,NR),lobatto_nodes(zs,NS),lobatto_nodes(zt,NT);
  array_init(struct pt_data,&testp,NEL*MULD(TN,TN,TN));
  crystal_init(&cr,&comm);
  crystal_init(&cr1,&comm1);
  test(&comm,&comm1);
  crystal_free(&cr1);
  crystal_free(&cr);
  array_free(&testp);
  
  comm_free(&comm);
  comm_free(&comm1);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}
