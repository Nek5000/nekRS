/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * The solver is based on Hiroshi Watanabe's 2D Gray-Scott reaction diffusion
 * code available at: https://github.com/kaityo256/sevendayshpc/tree/master/day5
 */

#include "gray-scott.h"

#include <mpi.h>
#include <stdexcept> // runtime_error
#include <vector>

GrayScott::GrayScott(const Settings &settings, MPI_Comm comm)
: settings(settings), comm(comm), rand_pool(5374857)
{
}

GrayScott::~GrayScott() {}

void GrayScott::init()
{
    init_mpi();
    init_field();
}

void GrayScott::iterate()
{
    auto temp_u = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, u);
    auto temp_v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, v);
    exchange(temp_u, temp_v);
    Kokkos::deep_copy(u, temp_u);
    Kokkos::deep_copy(v, temp_v);

    calc();

    std::swap(u, u2);
    std::swap(v, v2);
}

void GrayScott::restart(Kokkos::View<double ***, Kokkos::LayoutLeft> &u_in,
                        Kokkos::View<double ***, Kokkos::LayoutLeft> &v_in)
{
    auto const expected_len = (size_x + 2) * (size_y + 2) * (size_z + 2);
    if (u_in.size() == expected_len)
    {
        u = u_in;
        v = v_in;
    }
    else
    {
        throw std::runtime_error("Restart with incompatible array size, expected " +
                                 std::to_string(expected_len) + " got " +
                                 std::to_string(u_in.size()) + " elements");
    }
}

const Kokkos::View<double ***, Kokkos::LayoutLeft> GrayScott::u_ghost() const { return u; }

const Kokkos::View<double ***, Kokkos::LayoutLeft> GrayScott::v_ghost() const { return v; }

Kokkos::View<double ***, Kokkos::LayoutLeft> GrayScott::u_noghost() const
{
    return data_noghost(u);
}

Kokkos::View<double ***, Kokkos::LayoutLeft> GrayScott::v_noghost() const
{
    return data_noghost(v);
}

void GrayScott::u_noghost(Kokkos::View<double ***, Kokkos::LayoutLeft> u_no_ghost) const
{
    data_noghost(u, u_no_ghost);
}

void GrayScott::v_noghost(Kokkos::View<double ***, Kokkos::LayoutLeft> v_no_ghost) const
{
    data_noghost(v, v_no_ghost);
}

Kokkos::View<double ***, Kokkos::LayoutLeft>
GrayScott::data_noghost(const Kokkos::View<double ***, Kokkos::LayoutLeft> &data) const
{
    Kokkos::View<double ***, Kokkos::LayoutLeft> buf("noghost_temp", size_x, size_y, size_z);
    data_no_ghost_common(data, buf);
    return buf;
}

void GrayScott::data_noghost(const Kokkos::View<double ***, Kokkos::LayoutLeft> &data,
                             Kokkos::View<double ***, Kokkos::LayoutLeft> data_no_ghost) const
{
    data_no_ghost_common(data, data_no_ghost);
}

void GrayScott::init_field()
{
    Kokkos::resize(u, size_x + 2, size_y + 2, size_z + 2);
    Kokkos::deep_copy(u, 1.0);
    Kokkos::resize(v, size_x + 2, size_y + 2, size_z + 2);
    Kokkos::deep_copy(v, 0.0);
    Kokkos::resize(u2, size_x + 2, size_y + 2, size_z + 2);
    Kokkos::deep_copy(u2, 0.0);
    Kokkos::resize(v2, size_x + 2, size_y + 2, size_z + 2);
    Kokkos::deep_copy(v2, 0.0);

    const int d = 6;
    auto const L = settings.L;
    auto const settingsL = static_cast<int>(settings.L);
    auto const temp_u = u;
    auto const temp_v = v;
    size_t const ox = offset_x, oy = offset_y, oz = offset_z;
    size_t const sx = size_x, sy = size_y;
    auto const min_z = std::max(L / 2 - d, offset_z);
    auto const max_z = std::min(L / 2 + d, offset_z + size_z);
    Kokkos::parallel_for(
        "init_buffers", Kokkos::RangePolicy<>(min_z, max_z), KOKKOS_LAMBDA(int z) {
            for (int y = settingsL / 2 - d; y < settingsL / 2 + d; y++)
            {
                if (y < static_cast<int>(oy))
                    continue;
                if (y >= static_cast<int>(oy + sy))
                    continue;
                for (int x = settingsL / 2 - d; x < settingsL / 2 + d; x++)
                {
                    if (x < static_cast<int>(ox))
                        continue;
                    if (x >= static_cast<int>(ox + sx))
                        continue;
                    temp_u(x - ox + 1, y - oy + 1, z - oz + 1) = 0.25;
                    temp_v(x - ox + 1, y - oy + 1, z - oz + 1) = 0.33;
                }
            }
        });
}

void GrayScott::calc()
{
    auto const temp_u = u;
    auto const temp_v = v;
    auto const temp_u2 = u2;
    auto const temp_v2 = v2;
    auto const Du = settings.Du;
    auto const Dv = settings.Dv;
    auto const dt = settings.dt;
    auto const F = settings.F;
    auto const k = settings.k;
    auto const noise = settings.noise;
    size_t const sx = size_x, sy = size_y, sz = size_z;
    auto const random_pool = rand_pool;
    Kokkos::parallel_for(
        "calc_gray_scott", Kokkos::RangePolicy<>(1, sz + 1), KOKKOS_LAMBDA(int z) {
            RandomPool::generator_type generator = random_pool.get_state();
            double ts;
            for (int y = 1; y < static_cast<int>(sy) + 1; y++)
            {
                for (int x = 1; x < static_cast<int>(sx) + 1; x++)
                {
                    double du, dv;
                    // laplacian for u
                    ts = 0;
                    ts += temp_u(x - 1, y, z);
                    ts += temp_u(x + 1, y, z);
                    ts += temp_u(x, y - 1, z);
                    ts += temp_u(x, y + 1, z);
                    ts += temp_u(x, y, z - 1);
                    ts += temp_u(x, y, z + 1);
                    ts += -6.0 * temp_u(x, y, z);
                    ts /= 6.0;
                    du = Du * ts;

                    // laplacian for v
                    ts = 0;
                    ts += temp_v(x - 1, y, z);
                    ts += temp_v(x + 1, y, z);
                    ts += temp_v(x, y - 1, z);
                    ts += temp_v(x, y + 1, z);
                    ts += temp_v(x, y, z - 1);
                    ts += temp_v(x, y, z + 1);
                    ts += -6.0 * temp_v(x, y, z);
                    ts /= 6.0;
                    dv = Dv * ts;

                    du += (-temp_u(x, y, z) * temp_v(x, y, z) * temp_v(x, y, z) +
                           F * (1.0 - temp_u(x, y, z)));
                    dv += (temp_u(x, y, z) * temp_v(x, y, z) * temp_v(x, y, z) -
                           (F + k) * temp_v(x, y, z));
                    du += noise * generator.frand(-1.f, 1.f);
                    temp_u2(x, y, z) = temp_u(x, y, z) + du * dt;
                    temp_v2(x, y, z) = temp_v(x, y, z) + dv * dt;
                }
            }
            random_pool.free_state(generator);
        });
}

void GrayScott::init_mpi()
{
    int dims[3] = {};
    const int periods[3] = {1, 1, 1};
    int coords[3] = {};

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    MPI_Dims_create(procs, 3, dims);
    npx = dims[0];
    npy = dims[1];
    npz = dims[2];

    MPI_Cart_create(comm, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    px = coords[0];
    py = coords[1];
    pz = coords[2];

    size_x = settings.L / npx;
    size_y = settings.L / npy;
    size_z = settings.L / npz;

    if (px < settings.L % npx)
    {
        size_x++;
    }
    if (py < settings.L % npy)
    {
        size_y++;
    }
    if (pz < settings.L % npz)
    {
        size_z++;
    }

    offset_x = (settings.L / npx * px) + std::min(settings.L % npx, px);
    offset_y = (settings.L / npy * py) + std::min(settings.L % npy, py);
    offset_z = (settings.L / npz * pz) + std::min(settings.L % npz, pz);

    MPI_Cart_shift(cart_comm, 0, 1, &west, &east);
    MPI_Cart_shift(cart_comm, 1, 1, &down, &up);
    MPI_Cart_shift(cart_comm, 2, 1, &south, &north);

    // XY faces: size_x * (size_y + 2)
    MPI_Type_vector(static_cast<int>(size_y + 2), static_cast<int>(size_x),
                    static_cast<int>(size_x + 2), MPI_DOUBLE, &xy_face_type);
    MPI_Type_commit(&xy_face_type);

    // XZ faces: size_x * size_z
    MPI_Type_vector(static_cast<int>(size_z), static_cast<int>(size_x),
                    static_cast<int>((size_x + 2) * (size_y + 2)), MPI_DOUBLE, &xz_face_type);
    MPI_Type_commit(&xz_face_type);

    // YZ faces: (size_y + 2) * (size_z + 2)
    MPI_Type_vector(static_cast<int>((size_y + 2) * (size_z + 2)), 1, static_cast<int>(size_x + 2),
                    MPI_DOUBLE, &yz_face_type);
    MPI_Type_commit(&yz_face_type);
}

void GrayScott::exchange_xy(
    Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> local_data) const
{
    MPI_Status st;

    // Send XY face z=size_z to north and receive z=0 from south
    MPI_Sendrecv(&local_data.data()[l2i(1, 0, static_cast<int>(size_z))], 1, xy_face_type, north, 1,
                 &local_data.data()[l2i(1, 0, 0)], 1, xy_face_type, south, 1, cart_comm, &st);
    // Send XY face z=1 to south and receive z=size_z+1 from north
    MPI_Sendrecv(&local_data.data()[l2i(1, 0, 1)], 1, xy_face_type, south, 1,
                 &local_data.data()[l2i(1, 0, static_cast<int>(size_z + 1))], 1, xy_face_type,
                 north, 1, cart_comm, &st);
}

void GrayScott::exchange_xz(
    Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> local_data) const
{
    MPI_Status st;

    // Send XZ face y=size_y to up and receive y=0 from down
    MPI_Sendrecv(&local_data.data()[l2i(1, static_cast<int>(size_y), 1)], 1, xz_face_type, up, 2,
                 &local_data.data()[l2i(1, 0, 1)], 1, xz_face_type, down, 2, cart_comm, &st);
    // Send XZ face y=1 to down and receive y=size_y+1 from up
    MPI_Sendrecv(&local_data.data()[l2i(1, 1, 1)], 1, xz_face_type, down, 2,
                 &local_data.data()[l2i(1, static_cast<int>(size_y + 1), 1)], 1, xz_face_type, up,
                 2, cart_comm, &st);
}

void GrayScott::exchange_yz(
    Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> local_data) const
{
    MPI_Status st;

    // Send YZ face x=size_x to east and receive x=0 from west
    MPI_Sendrecv(&local_data.data()[l2i(static_cast<int>(size_x), 0, 0)], 1, yz_face_type, east, 3,
                 &local_data.data()[l2i(0, 0, 0)], 1, yz_face_type, west, 3, cart_comm, &st);
    // Send YZ face x=1 to west and receive x=size_x+1 from east
    MPI_Sendrecv(&local_data.data()[l2i(1, 0, 0)], 1, yz_face_type, west, 3,
                 &local_data.data()[l2i(static_cast<int>(size_x + 1), 0, 0)], 1, yz_face_type, east,
                 3, cart_comm, &st);
}

void GrayScott::exchange(Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> u,
                         Kokkos::View<double ***, Kokkos::LayoutLeft, Kokkos::HostSpace> v) const
{
    exchange_xy(u);
    exchange_xz(u);
    exchange_yz(u);

    exchange_xy(v);
    exchange_xz(v);
    exchange_yz(v);
}

void GrayScott::data_no_ghost_common(
    const Kokkos::View<double ***, Kokkos::LayoutLeft> &data,
    Kokkos::View<double ***, Kokkos::LayoutLeft> data_no_ghost) const
{
    auto const sx = size_x;
    auto const sy = size_y;
    auto const sz = size_z;
    Kokkos::parallel_for(
        "updateBuffer", Kokkos::RangePolicy<>(1, sz + 1), KOKKOS_LAMBDA(int z) {
            for (int y = 1; y < static_cast<int>(sy) + 1; y++)
            {
                for (int x = 1; x < static_cast<int>(sx) + 1; x++)
                {
                    data_no_ghost(x - 1, y - 1, z - 1) = data(x, y, z);
                }
            }
        });
}
