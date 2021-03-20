#include "genmap-impl.h"

int genmap_comm_size(genmap_comm c) { return (int)c->np; }

int genmap_comm_rank(genmap_comm c) { return (int)c->id; }

genmap_comm genmap_local_comm(genmap_handle h) { return h->local; }

genmap_comm genmap_global_comm(genmap_handle h) { return h->global; }

void genmap_comm_split(struct comm *old, int bin, int key, struct comm *new_) {
#ifdef MPI
  MPI_Comm new_comm;
  MPI_Comm_split(old->c, bin, key, &new_comm);
  comm_init(new_, new_comm);
  MPI_Comm_free(&new_comm);
#else
  comm_init(new_, 1);
#endif
}

void genmap_comm_scan(genmap_handle h, struct comm *c) {
  GenmapLong out[2][1], buf[2][1];
  GenmapLong lelt = genmap_get_nel(h);
  comm_scan(out, c, gs_long_long, gs_add, &lelt, 1, buf);
  genmap_set_local_start_index(h, out[0][0]);
  genmap_set_partition_nel(h, out[1][0]);
}
