#include <adios2.h>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <thread>
#include <vector>

#include <adios2/cxx11/KokkosView.h>

#include <Kokkos_Core.hpp>

int mpiRank, mpiSize;

template <class T, class MemSpace>
void PrintData(Kokkos::View<T *, MemSpace> &gpuData, const size_t step)
{
    auto data = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, gpuData);
    std::cout << "Rank: " << mpiRank << " Step: " << step << " [";
    for (int i = 0; i < data.extent_int(0); ++i)
    {
        std::cout << data(i) << " ";
    }
    std::cout << "]" << std::endl;
}

int main(int argc, char *argv[])
{
    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    // initialize adios2
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO dataManIO = adios.DeclareIO("whatever");
    dataManIO.SetEngine("DataMan");
    dataManIO.SetParameters({{"IPAddress", "127.0.0.1"}, {"Port", "12306"}, {"Timeout", "5"}});

    // open stream
    adios2::Engine dataManReader = dataManIO.Open("HelloDataMan", adios2::Mode::Read);

    // define variable
    adios2::Variable<float> floatArrayVar;

    Kokkos::DefaultExecutionSpace exe_space;
    std::cout << "Read on memory space: " << exe_space.name() << std::endl;
    // read data
    while (true)
    {
        auto status = dataManReader.BeginStep();
        if (status == adios2::StepStatus::OK)
        {
            floatArrayVar = dataManIO.InquireVariable<float>("FloatArray");
            auto shape = floatArrayVar.Shape();
            size_t datasize =
                std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
            Kokkos::View<float *, Kokkos::DefaultExecutionSpace::memory_space> floatVector(
                "simBuffer", datasize);
            dataManReader.Get<float>(floatArrayVar, floatVector, adios2::Mode::Sync);
            dataManReader.EndStep();
            PrintData(floatVector, dataManReader.CurrentStep());
        }
        else if (status == adios2::StepStatus::EndOfStream)
        {
            std::cout << "End of stream" << std::endl;
            break;
        }
    }

    // clean up
    dataManReader.Close();
    MPI_Finalize();

    return 0;
}
