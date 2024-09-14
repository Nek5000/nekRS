/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#ifndef __GRAY_SCOTT_H__
#define __GRAY_SCOTT_H__

#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <mpi.h>

#include "settings.h"

class GrayScott
{
public:
    // Dimension of process grid
    size_t npx, npy, npz;
    // Coordinate of this rank in process grid
    size_t px, py, pz;
    // Dimension of local array
    size_t size_x, size_y, size_z;
    // Offset of local array in the global array
    size_t offset_x, offset_y, offset_z;

    GrayScott(const Settings &settings, MPI_Comm comm);
    ~GrayScott();

    void init();
    void iterate();
    void restart(Kokkos::View<double ***, Kokkos::LayoutLeft> &u,
                 Kokkos::View<double ***, Kokkos::LayoutLeft> &v);

    const Kokkos::View<double ***, Kokkos::LayoutLeft> u_ghost() const;
    const Kokkos::View<double ***, Kokkos::LayoutLeft> v_ghost() const;

    Kokkos::View<double ***, Kokkos::LayoutLeft> u_noghost() const;
    Kokkos::View<double ***, Kokkos::LayoutLeft> v_noghost() const;

    void u_noghost(Kokkos::View<double ***, Kokkos::LayoutLeft> u_no_ghost) const;
    void v_noghost(Kokkos::View<double ***, Kokkos::LayoutLeft> v_no_ghost) const;

    Settings settings;

    Kokkos::View<double ***, Kokkos::LayoutLeft> u, v, u2, v2;

    int rank, procs;
    int west, east, up, down, north, south;
    MPI_Comm comm;
    MPI_Comm cart_comm;

    // MPI datatypes for halo exchange
    MPI_Datatype xy_face_type;
    MPI_Datatype xz_face_type;
    MPI_Datatype yz_face_type;

    using RandomPool = Kokkos::Random_XorShift64_Pool<Kokkos::DefaultExecutionSpace>;
    RandomPool rand_pool;

    // Setup cartesian communicator data types
    void init_mpi();
    // Setup initial conditions
    void init_field();

    // Process simulation for one timestep
    void calc();

    // Exchange faces with neighbors
    void exchange(Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> u,
                  Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> v) const;
    // Exchange XY faces with north/south
    void
    exchange_xy(Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> local_data) const;
    // Exchange XZ faces with up/down
    void
    exchange_xz(Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> local_data) const;
    // Exchange YZ faces with west/east
    void
    exchange_yz(Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> local_data) const;

    // Return a copy of data with ghosts removed
    Kokkos::View<double ***, Kokkos::LayoutLeft>
    data_noghost(const Kokkos::View<double ***, Kokkos::LayoutLeft> &data) const;

    // pointer version
    void data_noghost(const Kokkos::View<double ***, Kokkos::LayoutLeft> &data,
                      Kokkos::View<double ***, Kokkos::LayoutLeft> no_ghost) const;

    // Convert local coordinate to local index
    KOKKOS_FUNCTION int l2i(int x, int y, int z) const
    {
        return static_cast<int>(x + y * (size_x + 2) + z * (size_x + 2) * (size_y + 2));
    }

    void data_no_ghost_common(const Kokkos::View<double ***, Kokkos::LayoutLeft> &data,
                              Kokkos::View<double ***, Kokkos::LayoutLeft> data_no_ghost) const;
};

#endif
