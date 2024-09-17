/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <fstream>
#include <iostream>

#include <adios2.h>
#include <mpi.h>

#include "../common/timer.hpp"
#include "gray-scott.h"
#include "restart.h"
#include "writer.h"

void print_io_settings(const adios2::IO &io)
{
    std::cout << "Simulation writes data using engine type:              " << io.EngineType()
              << std::endl;
    auto ioparams = io.Parameters();
    std::cout << "IO parameters:  " << std::endl;
    for (const auto &p : ioparams)
    {
        std::cout << "    " << p.first << " = " << p.second << std::endl;
    }
}

void print_settings(const Settings &s, int restart_step)
{
    std::cout << "grid:             " << s.L << "x" << s.L << "x" << s.L << std::endl;
    if (restart_step > 0)
    {
        std::cout << "restart:          from step " << restart_step << std::endl;
    }
    else
    {
        std::cout << "restart:          no" << std::endl;
    }
    std::cout << "steps:            " << s.steps << std::endl;
    std::cout << "plotgap:          " << s.plotgap << std::endl;
    std::cout << "F:                " << s.F << std::endl;
    std::cout << "k:                " << s.k << std::endl;
    std::cout << "dt:               " << s.dt << std::endl;
    std::cout << "Du:               " << s.Du << std::endl;
    std::cout << "Dv:               " << s.Dv << std::endl;
    std::cout << "noise:            " << s.noise << std::endl;
    std::cout << "output:           " << s.output << std::endl;
    std::cout << "adios_config:     " << s.adios_config << std::endl;
}

void print_simulator_settings(const GrayScott &s)
{
    std::cout << "process layout:   " << s.npx << "x" << s.npy << "x" << s.npz << std::endl;
    std::cout << "local grid size:  " << s.size_x << "x" << s.size_y << "x" << s.size_z
              << std::endl;
}

int main(int argc, char **argv)
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 1;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    if (argc < 2)
    {
        if (rank == 0)
        {
            std::cerr << "Too few arguments" << std::endl;
            std::cerr << "Usage: gray-scott settings.json" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    Settings settings = Settings::from_json(argv[1]);

    GrayScott sim(settings, comm);
    sim.init();

    adios2::ADIOS adios(settings.adios_config, comm);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");
    adios2::IO io_ckpt = adios.DeclareIO("SimulationCheckpoint");

    int restart_step = 0;
    if (settings.restart)
    {
        restart_step = ReadRestart(comm, settings, sim, io_ckpt);
        io_main.SetParameter("AppendAfterSteps", std::to_string(restart_step / settings.plotgap));
    }

    Writer writer_main(settings, sim, io_main);
    writer_main.open(settings.output, (restart_step > 0));

    if (rank == 0)
    {
        print_io_settings(io_main);
        std::cout << "========================================" << std::endl;
        print_settings(settings, restart_step);
        print_simulator_settings(sim);
        std::cout << "========================================" << std::endl;
    }

#ifdef ENABLE_TIMERS
    Timer timer_total;
    Timer timer_compute;
    Timer timer_write;

    std::ostringstream log_fname;
    log_fname << "gray_scott_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
    log << "step\ttotal_gs\tcompute_gs\twrite_gs" << std::endl;
#endif

    for (int it = restart_step; it < settings.steps;)
    {
#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_total.start();
        timer_compute.start();
#endif

        sim.iterate();
        it++;

#ifdef ENABLE_TIMERS
        timer_compute.stop();
        MPI_Barrier(comm);
        timer_write.start();
#endif

        if (it % settings.plotgap == 0)
        {
            if (rank == 0)
            {
                std::cout << "Simulation at step " << it << " writing output step     "
                          << it / settings.plotgap << std::endl;
            }

            writer_main.write(it, sim);
        }

        if (settings.checkpoint && (it % settings.checkpoint_freq) == 0)
        {
            WriteCkpt(comm, it, settings, sim, io_ckpt);
        }

#ifdef ENABLE_TIMERS
        double time_write = timer_write.stop();
        double time_step = timer_total.stop();
        MPI_Barrier(comm);

        log << it << "\t" << timer_total.elapsed() << "\t" << timer_compute.elapsed() << "\t"
            << timer_write.elapsed() << std::endl;
#endif
    }

    writer_main.close();

#ifdef ENABLE_TIMERS
    log << "total\t" << timer_total.elapsed() << "\t" << timer_compute.elapsed() << "\t"
        << timer_write.elapsed() << std::endl;

    log.close();
#endif

    MPI_Finalize();
}
