#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

int log2i(sint i){
  sint k=1,l=0;
  while(k<=i) k*=2,l++;
  return l-1;
}

void setOwner(char *ptr,sint n,size_t inOffset,size_t outOffset,
  slong lelg,sint np)
{
  sint lelt=lelg/np;
  sint nrem=lelg%np;

  ulong *inPtr;
  sint  *outPtr;
  sint i; slong row;
  for(i=0; i<n; i++){
    inPtr =(ulong*)GETPTR(ptr,i,inOffset );
    outPtr=(sint  *)GETPTR(ptr,i,outOffset);
    row   =*inPtr-1;
    //FIXME: Assumes the 'reverse-Nek' element distribution
#if 0
    if(row<lelt*(np-nrem)) *outPtr=(sint) row/lelt;
    else *outPtr=np-nrem+(sint) (row-lelt*(np-nrem))/(lelt+1);
#else
    if(nrem==0) *outPtr=(sint) row/lelt;
    else if(row<(lelt+1)*nrem) *outPtr=(sint) row/(lelt+1);
    else *outPtr=nrem+(sint) (row-(lelt+1)*nrem)/lelt;
#endif
  }
}

// Following two functions can be combined
void compress_col(struct array *entries){
  GenmapScalar v; sint i,j;

  i=0; entry *ptr=entries->ptr; while(i<entries->n){
    v=ptr[i].v,j=i+1;
    while(j<entries->n && ptr[j].r==ptr[i].r && ptr[j].cn==ptr[i].cn)
      v+=ptr[j].v,ptr[j].c=0,j++;
    ptr[i].v=v;
    i=j;
  }
}

void mgLevelSetup(mgData d,uint lvl)
{
  assert(lvl>0); csr_mat M0=d->levels[lvl-1]->M;
  uint rn0=M0->rn,nnz0=M0->row_off[rn0];

  struct array entries=null_array;
  array_init(entry,&entries,nnz0);
  entries.n=nnz0;

  uint i,j,nn=0;
  entry *ptr=entries.ptr;
  for(i=0; i<rn0; i++)
    for(j=M0->row_off[i]; j<M0->row_off[i+1]; j++){
      ptr[nn].r=M0->row_start+i;
      ptr[nn].c=M0->col[j];
      ptr[nn].rn=(ptr[nn].r+1)/2;
      ptr[nn].cn=(ptr[nn].c+1)/2; // Let's collapse columns first
      ptr[nn].v=M0->v[j];
      nn++;
    }
  assert(nn==nnz0);

  slong out[2][1],bf[2][1],in=rn0;
  comm_scan(out,&d->c,gs_long,gs_add,&in,1,bf);
  slong ng=out[1][0];

  ulong ngc=ng/2;
  if(ngc==0) return;

  if(ng>1 && ng%2==1)
    for(j=0; j<M0->row_off[rn0]; j++){
      if(ptr[j].c==ng) ptr[j].cn-=1;
      if(ptr[j].r==ng) ptr[j].rn-=1;
    }

  uint npc=d->c.np;
  if(ngc<npc) npc=ngc;

  /* setup gs ids for fine level (rhs interpolation) */
  ptr=entries.ptr;
  slong *ids; GenmapMalloc(rn0,&ids);
  for(i=j=0; i<nnz0; i++)
    if(ptr[i].r==ptr[i].c) ids[j++]=-ptr[i].cn;
  assert(j==rn0);

  /* coarsen the cols */
  buffer buf; buffer_init(&buf,1024);
  if(entries.n){
    sarray_sort_2(entry,entries.ptr,entries.n,r,1,cn,1,&buf);
    compress_col(&entries);
  }

  struct crystal cr;
  crystal_init(&cr,&d->c);

  setOwner(entries.ptr,nnz0,offsetof(entry,rn),offsetof(entry,p),ngc,npc);
  sarray_transfer(entry,&entries,p,1,&cr);

  // sort by rn and cn
  sarray_sort_2(entry,entries.ptr,entries.n,rn,1,cn,1,&buf);

  i=j=nn=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn) j++;
    i=j,nn++;
  }

  /* create the matrix */
  GenmapMalloc(1,&d->levels[lvl]); mgLevel l=d->levels[lvl];l->data=d;

  GenmapMalloc(1,&l->M); csr_mat M1 =l->M; M1->rn=nn;
  GenmapMalloc(M1->rn+1,&M1->row_off);

  slong cn=nn; comm_scan(out,&d->c,gs_long,gs_add,&cn,1,bf);
  M1->row_start=out[0][0]+1;

  uint nnz1=0;
  i=j=0; ptr=entries.ptr; while(i<entries.n){
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn)
      j++;
    i=j,nnz1++;
  }

  if(nnz1==0)
    M1->col=NULL,M1->v=NULL,M1->diag=NULL;
  else{
    GenmapMalloc(nnz1,&M1->col),GenmapMalloc(nnz1,&M1->v);
    GenmapMalloc(M1->rn,&M1->diag);
  }

  uint rn1; GenmapScalar v;
  M1->row_off[0]=i=j=nn=rn1=0; ptr=entries.ptr; while(i<entries.n){
    v=0.0;
    while(j<entries.n && ptr[j].rn==ptr[i].rn && ptr[j].cn==ptr[i].cn){
      if(ptr[j].c>0) v+=ptr[j].v;
      j++;
    }
    M1->col[nn]=ptr[i].cn,M1->v[nn]=v,nn++;

    if((j<entries.n && ptr[j].rn!=ptr[i].rn) || j>=entries.n)
      M1->row_off[++rn1]=nn;
    i=j;
  }
  assert(nn==nnz1); //sanity check
  assert(rn1==M1->rn); //sanity check

  /* setup gs ids for coarse level (rhs interpolation ) */
  GenmapRealloc(rn0+rn1,&ids);
  for(i=nn=0; i<rn1; i++)
    for(j=M1->row_off[i]; j<M1->row_off[i+1]; j++)
      if(M1->row_start+i==M1->col[j]){
        ids[rn0+nn]=M1->col[j];
        M1->diag[i]=M1->v[j];
        nn++;
      }
  assert(nn==M1->rn);

  d->levels[lvl-1]->J=gs_setup(ids,rn0+M1->rn,&d->c,0,gs_crystal_router,0);

  /* setup gs handle for the mat-vec */
  GenmapRealloc(nnz1,&ids);
  for(i=0; i<M1->rn; i++)
    for(j=M1->row_off[i]; j<M1->row_off[i+1]; j++)
      if(M1->row_start+i==M1->col[j]) ids[j]=M1->col[j];
      else ids[j]=-M1->col[j];

  M1->gsh=gs_setup(ids,nnz1,&d->c,0,gs_crystal_router,0);

  GenmapFree(ids);
  buffer_free(&buf);
  crystal_free(&cr);
  array_free(&entries);
}

void mgSetup(GenmapComm c,csr_mat M,mgData *d_){
  GenmapMalloc(1,d_); mgData d=*d_; comm_dup(&d->c,&c->gsc);

  uint np=GenmapCommSize(c); uint rn=M->rn;

  slong out[2][1],bf[2][1],in=rn;
  comm_scan(out,&d->c,gs_long,gs_add,&in,1,bf);
  slong rg=out[1][0];

  d->nlevels=log2i(rg)+1;
  GenmapMalloc(d->nlevels  ,&d->levels   );
  GenmapMalloc(d->nlevels+1,&d->level_off);

  GenmapMalloc(1,&d->levels[0]);
  d->levels[0]->M=M; d->level_off[0]=0; d->level_off[1]=M->rn;
  d->levels[0]->nsmooth=2,d->levels[0]->sigma=0.6;

  uint i; uint nnz=M->row_off[M->rn];
  for(i=1;i<d->nlevels;i++){
    mgLevelSetup(d,i); csr_mat Mi=d->levels[i]->M;
    if(Mi->row_off[Mi->rn]>nnz)
      nnz=Mi->row_off[Mi->rn];
    d->level_off[i+1]=d->level_off[i]+Mi->rn;
    d->levels[i]->nsmooth=2;
    d->levels[i]->sigma  =0.6;
  }

  GenmapMalloc(d->level_off[d->nlevels],&d->x  );
  GenmapMalloc(d->level_off[d->nlevels],&d->y  );
  GenmapMalloc(d->level_off[d->nlevels],&d->b  );
  GenmapMalloc(d->level_off[d->nlevels],&d->u  );
  GenmapMalloc(d->level_off[d->nlevels],&d->rhs);
  GenmapMalloc(nnz                     ,&d->buf);
}

void mgFree(mgData d){
  mgLevel *l=d->levels;
  uint i,nlevels=d->nlevels;
  for(i=0; i<nlevels; i++){
    if(i>0)
      csr_mat_free(l[i]->M);
    if(i<nlevels-1){
      gs_free(l[i]->J); GenmapFree(l[i]);
    }
  }

  GenmapFree(l);
  GenmapFree(d->level_off);
  GenmapFree(d->y); GenmapFree(d->x); GenmapFree(d->b);
  GenmapFree(d->buf); GenmapFree(d->rhs); GenmapFree(d->u);
}
