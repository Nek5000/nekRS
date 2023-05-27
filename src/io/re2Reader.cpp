#include "nrs.hpp"
#include "re2Reader.hpp"

void re2::nelg(const std::string& meshFile, int& nelgt, int& nelgv, MPI_Comm comm)
{
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  if(rank == 0) {
    const int re2HeaderBytes = 80;
    char *buf = (char*) calloc(std::max(re2HeaderBytes, (int)meshFile.length()+1), sizeof(char));
    strcpy(buf, meshFile.c_str());
    FILE *fp = fopen(buf, "r");
    nrsCheck(!fp, MPI_COMM_SELF, EXIT_FAILURE, "Cannot find %s!\n", buf);
    fgets(buf, re2HeaderBytes, fp);
    fclose(fp);
 
    char ver[6];
    sscanf(buf, "%5s", ver);

    int ndim;
    if(strcmp(ver, "#v004") == 0) {
      sscanf(buf, "%5s %d %d %d", ver, &nelgt, &ndim, &nelgv);
    } else if(strcmp(ver, "#v002") == 0 || strcmp(ver, "#v003") == 0) { 
      sscanf(buf, "%5s %9d %1d %9d", ver, &nelgt, &ndim, &nelgv); 
    } else {
      nrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "Unsupported re2 version %5s!\n", ver);
    }

    nrsCheck(ndim != 3, MPI_COMM_SELF, EXIT_FAILURE,
             "\nUnsupported ndim=%d read from re2 header!\n", ndim);
    
    nrsCheck(nelgt <= 0 || nelgv <=0 || nelgv > nelgt, MPI_COMM_SELF, EXIT_FAILURE,
             "\nInvalid nelgt=%d / nelgv=%d read from re2 header!\n", nelgt, nelgv);

    free(buf);
  }

  MPI_Bcast(&nelgt, 1, MPI_INT, 0, comm);
  MPI_Bcast(&nelgv, 1, MPI_INT, 0, comm);
}
