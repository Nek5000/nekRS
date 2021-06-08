
#define hash_data              TOKEN_PASTE(findpts_hash_data_               ,D)
#define eval_src_pt            TOKEN_PASTE(eval_src_pt_                     ,D)
#define eval_out_pt            TOKEN_PASTE(eval_out_pt_                     ,D)
#define findpts_el_data        TOKEN_PASTE(findpts_el_data_                 ,D)
#define findpts_el_eval        TOKEN_PASTE(ogs_findpts_el_eval_             ,D)
#define findpts_local_data     TOKEN_PASTE(findpts_local_data_              ,D)
#define findpts_local_eval     TOKEN_PASTE(ogs_findpts_local_eval_          ,D)
#define findpts_data           TOKEN_PASTE(findpts_data_                    ,D)
#define findpts_eval           TOKEN_PASTE(ogs_findpts_eval_                ,D)
#define cpp_findpts_local_eval TOKEN_PASTE(ogs_findpts_local_eval_internal_ ,D)

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

struct findpts_data {
  struct crystal cr;
  struct findpts_local_data local;
  struct hash_data hash;
};

struct eval_src_pt { double r[D]; uint index, proc, el; };
struct eval_out_pt { double out; uint index, proc; };

void findpts_local_eval(      double *const out_base, const unsigned out_stride,
                        const uint   *const  el_base, const unsigned  el_stride,
                        const double *const   r_base, const unsigned   r_stride,
                        const uint npt,
                        const void *const in, struct findpts_local_data *const fd)
{
  unsigned lag_data_size[D];
  for (int i = 0; i < D; ++i) lag_data_size[i] = gll_lag_size(fd->fed.n[i]);
  cpp_findpts_local_eval(out_base, out_stride,
                          el_base,  el_stride,
                           r_base,   r_stride,
                         npt, in, fd->ntot,
                         fd->fed.n, fd->fed.lag_data, lag_data_size);
}

void findpts_eval(      double *const  out_base, const unsigned  out_stride,
                  const uint   *const code_base, const unsigned code_stride,
                  const uint   *const proc_base, const unsigned proc_stride,
                  const uint   *const   el_base, const unsigned   el_stride,
                  const double *const    r_base, const unsigned    r_stride,
                  const uint npt,
                  const void *const in, struct findpts_data *const fd)
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
    findpts_local_eval(&opt->out ,sizeof(struct eval_out_pt),
                       &spt->el  ,sizeof(struct eval_src_pt),
                        spt->r   ,sizeof(struct eval_src_pt),
                       src.n, in,&fd->local);
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

#undef findpts_local_eval
#undef findpts_local_data
#undef findpts_el_eval
#undef findpts_el_data
