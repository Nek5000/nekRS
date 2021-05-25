/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef OGS_INTERFACE_H
#define OGS_INTERFACE_H 1

extern "C"
{
  void *ogsHostSetup(MPI_Comm comm, dlong Ngather, hlong *gatherIds, int unique, int verbose);
  void  ogsGsUnique(hlong *gatherIds, dlong Ngather, MPI_Comm comm);

  void ogsHostGatherScatter    (void *v, const char *type, const char *op, void *gsh);
  void ogsHostGatherScatterVec (void *v, const int k, const char *type, const char *op, void *gsh);
  void ogsHostGatherScatterMany(void *v, const int k, const char *type, const char *op, void *gsh);

  void ogsHostGather    (void *v, const char *type, const char *op, void *gsh);
  void ogsHostGatherVec (void *v, const int k, const char *type, const char *op, void *gsh);
  void ogsHostGatherMany(void *v, const int k, const char *type, const char *op, void *gsh);

  void ogsHostScatter    (void *v, const char *type, const char *op, void *gsh);
  void ogsHostScatterVec (void *v, const int k, const char *type, const char *op, void *gsh);
  void ogsHostScatterMany(void *v, const int k, const char *type, const char *op, void *gsh);

  struct findpts_data_2;
  struct findpts_data_3;

  typedef findpts_data_2 ogsFindptsData_2;
  typedef findpts_data_3 ogsFindptsData_3;

  struct findpts_data_2 *ogsHostFindptsSetup_2(
    MPI_Comm comm,
    const dfloat *const elx[2],
    const dlong n[2], const dlong nel,
    const dlong m[2], const dfloat bbox_tol,
    const hlong local_hash_size, const hlong global_hash_size,
    const dlong npt_max, const dfloat newt_tol);

  struct findpts_data_3 *ogsHostFindptsSetup_3(
    MPI_Comm comm,
    const dfloat *const elx[3],
    const dlong n[3], const dlong nel,
    const dlong m[3], const dfloat bbox_tol,
    const hlong local_hash_size, const hlong global_hash_size,
    const dlong npt_max, const dfloat newt_tol);


  void ogsHostFindptsFree_2(struct findpts_data_2 *fd);
  void ogsHostFindptsFree_3(struct findpts_data_3 *fd);

  void ogsHostFindpts_2(    dlong  *const  code_base   , const dlong  code_stride   ,
                            dlong  *const  proc_base   , const dlong  proc_stride   ,
                            dlong  *const    el_base   , const dlong    el_stride   ,
                            dfloat *const     r_base   , const dlong     r_stride   ,
                            dfloat *const dist2_base   , const dlong dist2_stride   ,
                      const dfloat *const     x_base[2], const dlong     x_stride[2],
                      const dfloat npt, struct findpts_data_2 *const fd);

  void ogsHostFindpts_3(    dlong  *const  code_base   , const dlong  code_stride   ,
                            dlong  *const  proc_base   , const dlong  proc_stride   ,
                            dlong  *const    el_base   , const dlong    el_stride   ,
                            dfloat *const     r_base   , const dlong     r_stride   ,
                            dfloat *const dist2_base   , const dlong dist2_stride   ,
                      const dfloat *const     x_base[3], const dlong     x_stride[3],
                      const dfloat npt, struct findpts_data_3 *const fd);

  void ogsHostFindptsEval_2(
          dfloat *const  out_base, const dlong  out_stride,
    const dlong  *const code_base, const dlong code_stride,
    const dlong  *const proc_base, const dlong proc_stride,
    const dlong  *const   el_base, const dlong   el_stride,
    const dfloat *const    r_base, const dlong    r_stride,
    const dlong npt, const dfloat *const in, struct findpts_data_2 *const fd);

  void ogsHostFindptsEval_3(
          dfloat *const  out_base, const dlong  out_stride,
    const dlong  *const code_base, const dlong code_stride,
    const dlong  *const proc_base, const dlong proc_stride,
    const dlong  *const   el_base, const dlong   el_stride,
    const dfloat *const    r_base, const dlong    r_stride,
    const dlong npt, const dfloat *const in, struct findpts_data_3 *const fd);

  void ogsDevFindptsEval_2(
          dfloat *const  out_base, const dlong  out_stride,
    const dlong  *const code_base, const dlong code_stride,
    const dlong  *const proc_base, const dlong proc_stride,
    const dlong  *const   el_base, const dlong   el_stride,
    const dfloat *const    r_base, const dlong    r_stride,
    const dlong npt, void *const in, struct findpts_data_2 *const fd);

  void ogsDevFindptsEval_3(
          dfloat *const  out_base, const dlong  out_stride,
    const dlong  *const code_base, const dlong code_stride,
    const dlong  *const proc_base, const dlong proc_stride,
    const dlong  *const   el_base, const dlong   el_stride,
    const dfloat *const    r_base, const dlong    r_stride,
    const dlong npt, void *const in, struct findpts_data_3 *const fd);

  void ogsHostFree(void *gsh);

}

#endif
