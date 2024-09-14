#include <iostream>
#include <vector>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
    int rank = 0;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    adios2::ADIOS ad;
    try
    {
#if ADIOS2_USE_MPI
        ad = adios2::ADIOS("does_not_exist.xml", MPI_COMM_WORLD);
#else
        ad = adios2::ADIOS("does_not_exist.xml");
#endif
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cout << e.what() << "\n";
        }
#if ADIOS2_USE_MPI
        ad = adios2::ADIOS(MPI_COMM_WORLD);
#else
        ad = adios2::ADIOS();
#endif
    }

#if ADIOS2_USE_MPI
    return MPI_Finalize();
#endif
}
