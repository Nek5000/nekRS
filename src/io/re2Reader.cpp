#include "nrs.hpp"
#include "re2Reader.hpp"

void re2::nelg(const std::string& meshFile, int& nelgt, int& nelgv, MPI_Comm comm)
{
  int rank = 0;
  if(comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);

  int err = 0;
  if(rank == 0) {
    char buf[FILENAME_MAX];
    strcpy(buf, meshFile.c_str());
    FILE *fp = fopen(buf, "r");
    if (!fp) {
      if(rank == 0) printf("\nERROR: Cannot find %s!\n", buf);
      ABORT(EXIT_FAILURE);
    }
    fgets(buf, 80, fp);
    fclose(fp);
 
    char ver[6];
    int ndim;
    // has to match header in re2
    sscanf(buf, "%5s %9d %1d %9d", ver, &nelgt, &ndim, &nelgv);

    if(ndim != 3) {
      if(rank == 0) printf("\nERROR: Unsupported ndim=%d read from re2 header!\n", ndim);
      err++;
    }
    if(nelgt <= 0 || nelgv <=0 || nelgv > nelgt) {
      if(rank == 0) printf("\nERROR: Invalid nelgt=%d / nelgv=%d read from re2 header!\n", nelgt, nelgv);
      err++;
    }
  }
  if(comm != MPI_COMM_NULL) MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, comm); 
  if(err) ABORT(EXIT_FAILURE);

  if(comm != MPI_COMM_NULL) MPI_Bcast(&nelgt, 1, MPI_INT, 0, comm);
  if(comm != MPI_COMM_NULL) MPI_Bcast(&nelgv, 1, MPI_INT, 0, comm);
}
