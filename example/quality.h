void printPartStat(long long *vtx, int nel, int nv, comm_ext ce) {
  int i, j;

  struct comm comm;
  int np, id;

  int Nmsg;
  int *Ncomm;

  int nelMin, nelMax;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax, nssSum;

  struct gs_data *gsh;
  int b;

  int numPoints;
  long long *data;

  comm_init(&comm, ce);
  np = comm.np;
  id = comm.id;

  if(np == 1) return;

  numPoints = nel * nv;
  data = (long long*) malloc(numPoints * sizeof(long long));
  for(i = 0; i < numPoints; i++) data[i] = vtx[i];

  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  pw_data_nmsg(gsh, &Nmsg);
  Ncomm = (int *) malloc(Nmsg * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for(i = 1; i < Nmsg; ++i) {
    nsMax = Ncomm[i] > Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsMin = Ncomm[i] < Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsSum += Ncomm[i];
  }
  comm_allreduce(&comm, gs_int, gs_max, &nsMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nsMin, 1, &b);

  nssMin = nsSum;
  nssMax = nsSum;
  nssSum = nsSum;
  comm_allreduce(&comm, gs_int, gs_max, &nssMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nssMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nssSum, 1, &b);

  nsSum = nsSum / Nmsg;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum, 1, &b);

  nelMax = nel;
  nelMin = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);

  if(id == 0) {
    printf(
      " Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
      ncMax, ncMin, (double)ncSum / np);
    printf(
      " Max nvolume: %d | Min nvolume: %d | Avg nvolume: %lf\n",
      nsMax, nsMin, (double)nsSum / np);
    printf(
      " Max volume: %d | Min volume: %d | Avg volume: %lf\n",
      nssMax, nssMin, (double)nssSum / np);
    printf(
      " Max elements: %d | Min elements: %d | Balance: %lf\n",
      nelMax, nelMin, (double)nelMax / nelMin);
    fflush(stdout);
  }

  comm_free(&comm);
}
