
#define hash_data              TOKEN_PASTE(findpts_hash_data_               ,D)
#define hash_index             TOKEN_PASTE(hash_index_                      ,D)
#define eval_src_pt            TOKEN_PASTE(eval_src_pt_                     ,D)
#define eval_out_pt            TOKEN_PASTE(eval_out_pt_                     ,D)
#define src_pt                 TOKEN_PASTE(src_pt_                          ,D)
#define out_pt                 TOKEN_PASTE(out_pt_                          ,D)
#define findpts_el_data        TOKEN_PASTE(findpts_el_data_                 ,D)
#define findpts_local_data     TOKEN_PASTE(findpts_local_data_              ,D)
#define findpts_local_eval_int TOKEN_PASTE(ogs_findpts_local_eval_internal_ ,D)
#define findpts_local          TOKEN_PASTE(ogs_findpts_local_               ,D)
#define findpts_local_eval     TOKEN_PASTE(ogs_findpts_local_eval_          ,D)
#define findpts_data           TOKEN_PASTE(findpts_data_                    ,D)
#define findpts                TOKEN_PASTE(ogs_findpts_                     ,D)
#define findpts_eval           TOKEN_PASTE(ogs_findpts_eval_                ,D)

#define   AT(T,var,i)   \
        (T*)(      (char*)var##_base   +(i)*var##_stride   )
#define  CAT(T,var,i) \
  (const T*)((const char*)var##_base   +(i)*var##_stride   )
#define CATD(T,var,i,d) \
  (const T*)((const char*)var##_base[d]+(i)*var##_stride[d])

struct hash_data {
  ulong hash_n;
  struct dbl_range bnd[D];
  double fac[D];
  uint *offset;
};

static ulong hash_index(const struct hash_data *p, const double x[D])
{
  const ulong n = p->hash_n;
  return ( WHEN_3D( hash_index_aux(p->bnd[2].min,p->fac[2],n,x[2])  *n )
                   +hash_index_aux(p->bnd[1].min,p->fac[1],n,x[1]) )*n
                   +hash_index_aux(p->bnd[0].min,p->fac[0],n,x[0]);
}


struct findpts_data {
  struct crystal cr;
  struct findpts_local_data local;
  struct hash_data hash;
};

struct src_pt { double x[D]; uint index, proc; };
struct out_pt { double r[D], dist2; uint index, code, el, proc; };

void findpts(      uint   *const  code_base   , const unsigned  code_stride   ,
                   uint   *const  proc_base   , const unsigned  proc_stride   ,
                   uint   *const    el_base   , const unsigned    el_stride   ,
                   double *const     r_base   , const unsigned     r_stride   ,
                   double *const dist2_base   , const unsigned dist2_stride   ,
             const double *const     x_base[D], const unsigned     x_stride[D],
             const uint npt, struct findpts_data *const fd,
             const void *const ogs_fd)
{
  const uint np = fd->cr.comm.np, id=fd->cr.comm.id;
  struct array hash_pt, src_pt, out_pt;
  /* look locally first */
  if(npt) findpts_local( code_base, code_stride,
                           el_base,   el_stride,
                            r_base,    r_stride,
                        dist2_base,dist2_stride,
                            x_base,    x_stride,
                        npt, ogs_fd);
  /* send unfound and border points to global hash cells */
  {
    uint index;
    uint *code=code_base, *proc=proc_base;
    const double *xp[D];
    struct src_pt *pt;
    unsigned d; for(d=0;d<D;++d) xp[d]=x_base[d];
    array_init(struct src_pt, &hash_pt, npt), pt=hash_pt.ptr;
    for(index=0;index<npt;++index) {
      double x[D]; for(d=0;d<D;++d) x[d]=*xp[d];
      *proc = id;
      if(*code!=CODE_INTERNAL) {
        const uint hi = hash_index(&fd->hash,x);
        unsigned d;
        for(d=0;d<D;++d) pt->x[d]=x[d];
        pt->index=index;
        pt->proc=hi%np;
        ++pt;
      }
      for(d=0;d<D;++d)
      xp[d] = (const double*)((const char*)xp[d]+   x_stride[d]);
      code  =         (uint*)(      (char*)code +code_stride   );
      proc  =         (uint*)(      (char*)proc +proc_stride   );
    }
    hash_pt.n = pt - (struct src_pt*)hash_pt.ptr;
    sarray_transfer(struct src_pt,&hash_pt,proc,1,&fd->cr);
  }
  /* look up points in hash cells, route to possible procs */
  {
    const uint *const hash_offset = fd->hash.offset;
    uint count=0, *proc, *proc_p;
    const struct src_pt *p = hash_pt.ptr, *const pe = p+hash_pt.n;
    struct src_pt *q;
    for(;p!=pe;++p) {
      const uint hi = hash_index(&fd->hash,p->x)/np;
      const uint i = hash_offset[hi], ie = hash_offset[hi+1];
      count += ie-i;
    }
    proc_p = proc = tmalloc(uint,count);
    array_init(struct src_pt,&src_pt,count), q=src_pt.ptr;
    for(p=hash_pt.ptr;p!=pe;++p) {
      const uint hi = hash_index(&fd->hash,p->x)/np;
      uint i = hash_offset[hi]; const uint ie = hash_offset[hi+1];
      for(;i!=ie;++i) {
        const uint pp = hash_offset[i];
        if(pp==p->proc) continue; /* don't send back to source proc */
        *proc_p++ = pp;
        *q++ = *p;
      }
    }
    array_free(&hash_pt);
    src_pt.n = proc_p-proc;
    #ifdef DIAGNOSTICS
    printf("(proc %u) hashed; routing %u/%u\n",id,(unsigned)src_pt.n,count);
    #endif
    sarray_transfer_ext(struct src_pt,&src_pt,proc,sizeof(uint),&fd->cr);
    free(proc);
  }
  /* look for other procs' points, send back */
  {
    uint n=src_pt.n;
    const struct src_pt *spt;
    struct out_pt *opt;
    array_init(struct out_pt,&out_pt,n), out_pt.n=n;
    spt=src_pt.ptr, opt=out_pt.ptr;
    for(;n;--n,++spt,++opt) opt->index=spt->index,opt->proc=spt->proc;
    spt=src_pt.ptr, opt=out_pt.ptr;
    if(src_pt.n) {
      const double *spt_x_base[D]; unsigned spt_x_stride[D];
      unsigned d; for(d=0;d<D;++d) spt_x_base[d] = spt[0].x+d,
                                   spt_x_stride[d] = sizeof(struct src_pt);
      findpts_local(&opt[0].code ,sizeof(struct out_pt),
                    &opt[0].el   ,sizeof(struct out_pt),
                     opt[0].r    ,sizeof(struct out_pt),
                    &opt[0].dist2,sizeof(struct out_pt),
                     spt_x_base  ,spt_x_stride,
                    src_pt.n, ogs_fd);
    }
    array_free(&src_pt);
    /* group by code to eliminate unfound points */
    sarray_sort(struct out_pt,opt,out_pt.n, code,0, &fd->cr.data);
    n=out_pt.n; while(n && opt[n-1].code==CODE_NOT_FOUND) --n;
    out_pt.n=n;
    #ifdef DIAGNOSTICS
    printf("(proc %u) sending back %u found points\n",id,(unsigned)out_pt.n);
    #endif
    sarray_transfer(struct out_pt,&out_pt,proc,1,&fd->cr);
  }
  /* merge remote results with user data */
  {
    uint n=out_pt.n;
    struct out_pt *opt;
    for(opt=out_pt.ptr;n;--n,++opt) {
      const uint index = opt->index;
      uint *code = AT(uint,code,index);
      double *dist2 = AT(double,dist2,index);
      if(*code==CODE_INTERNAL) continue;
      if(*code==CODE_NOT_FOUND
         || opt->code==CODE_INTERNAL
         || opt->dist2<*dist2) {
        double *r = AT(double,r,index);
        uint  *el = AT(uint,el,index), *proc = AT(uint,proc,index);
        unsigned d; for(d=0;d<D;++d) r[d]=opt->r[d];
        *dist2 = opt->dist2;
        *proc = opt->proc;
        *el = opt->el;
        *code = opt->code;
      }
    }
    array_free(&out_pt);
  }
}


void findpts_eval(      double *const  out_base, const unsigned  out_stride,
                  const uint   *const code_base, const unsigned code_stride,
                  const uint   *const proc_base, const unsigned proc_stride,
                  const uint   *const   el_base, const unsigned   el_stride,
                  const double *const    r_base, const unsigned    r_stride,
                  const uint npt,
                  const void *const in, struct findpts_data *const fd,
                  const void *const ogs_fd)
{
  struct array src, outpt;
  /* copy user data, weed out unfound points, send out */
  {
    uint index;
    const uint *code=code_base, *proc=proc_base, *el=el_base;
    const double *r=r_base;
    struct eval_src_pt *pt;
    array_init(struct eval_src_pt, &src, npt), pt=src.ptr;
    for(index=0;index<npt;++index) {
      if(*code!=CODE_NOT_FOUND) {
        unsigned d;
        for(d=0;d<D;++d) pt->r[d]=r[d];
        pt->index=index;
        pt->proc=*proc;
        pt->el=*el;
        ++pt;
      }
      r    = (const double*)((const char*)r   +   r_stride);
      code = (const   uint*)((const char*)code+code_stride);
      proc = (const   uint*)((const char*)proc+proc_stride);
      el   = (const   uint*)((const char*)el  +  el_stride);
    }
    src.n = pt - (struct eval_src_pt*)src.ptr;
    sarray_transfer(struct eval_src_pt,&src,proc,1,&fd->cr);
  }
  /* evaluate points, send back */
  {
    uint n=src.n;
    const struct eval_src_pt *spt;
    struct eval_out_pt *opt;
    /* group points by element */
    sarray_sort(struct eval_src_pt,src.ptr,n, el,0, &fd->cr.data);
    array_init(struct eval_out_pt,&outpt,n), outpt.n=n;
    spt=src.ptr, opt=outpt.ptr;
    findpts_local_eval_int(opt, spt, src.n, in,&fd->local, ogs_fd);
    spt=src.ptr, opt=outpt.ptr;
    for(;n;--n,++spt,++opt) opt->index=spt->index,opt->proc=spt->proc;
    array_free(&src);
    sarray_transfer(struct eval_out_pt,&outpt,proc,1,&fd->cr);
  }
  /* copy results to user data */
  {
    uint n=outpt.n;
    struct eval_out_pt *opt;
    for(opt=outpt.ptr;n;--n,++opt) *AT(double,out,opt->index)=opt->out;
    array_free(&outpt);
  }
}

#undef CATD
#undef CAT
#undef AT

#undef hash_data
#undef eval_src_pt
#undef eval_out_pt
#undef src_pt
#undef out_pt
#undef findpts_el_data
#undef findpts_local_data
#undef findpts_local
#undef findpts_local_eval
#undef findpts_data
#undef findpts
#undef findpts_eval
