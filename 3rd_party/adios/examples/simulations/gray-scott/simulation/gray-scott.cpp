/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * The solver is based on Hiroshi Watanabe's 2D Gray-Scott reaction diffusion
 * code available at:
 * https://github.com/kaityo256/sevendayshpc/tree/master/day5
 */

#include "gray-scott.h"

#include <mpi.h>
#include <random>
#include <stdexcept> // runtime_error
#include <vector>

GrayScott::GrayScott(const Settings &settings, MPI_Comm comm)
: settings(settings), comm(comm), rand_dev(), mt_gen(rand_dev()), uniform_dist(-1.0, 1.0)
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
    exchange(u, v);
    calc(u, v, u2, v2);

    u.swap(u2);
    v.swap(v2);
}

void GrayScott::restart(std::vector<double> &u_in, std::vector<double> &v_in)
{
    auto expected_len = (size_x + 2) * (size_y + 2) * (size_z + 2);
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

const std::vector<double> &GrayScott::u_ghost() const { return u; }

const std::vector<double> &GrayScott::v_ghost() const { return v; }

std::vector<double> GrayScott::u_noghost() const { return data_noghost(u); }

std::vector<double> GrayScott::v_noghost() const { return data_noghost(v); }

void GrayScott::u_noghost(double *u_no_ghost) const { data_noghost(u, u_no_ghost); }

void GrayScott::v_noghost(double *v_no_ghost) const { data_noghost(v, v_no_ghost); }

std::vector<double> GrayScott::data_noghost(const std::vector<double> &data) const
{
    std::vector<double> buf(size_x * size_y * size_z);
    data_no_ghost_common(data, buf.data());
    return buf;
}

void GrayScott::data_noghost(const std::vector<double> &data, double *data_no_ghost) const
{
    data_no_ghost_common(data, data_no_ghost);
}

void GrayScott::init_field()
{
    const int V = static_cast<int>((size_x + 2) * (size_y + 2) * (size_z + 2));
    u.resize(V, 1.0);
    v.resize(V, 0.0);
    u2.resize(V, 0.0);
    v2.resize(V, 0.0);

    const int d = 6;
    const auto settingsL = static_cast<int>(settings.L);
    for (int z = settingsL / 2 - d; z < settingsL / 2 + d; z++)
    {
        for (int y = settingsL / 2 - d; y < settingsL / 2 + d; y++)
        {
            for (int x = settingsL / 2 - d; x < settingsL / 2 + d; x++)
            {
                if (!is_inside(x, y, z))
                    continue;
                int i = g2i(x, y, z);
                u[i] = 0.25;
                v[i] = 0.33;
            }
        }
    }
}

double GrayScott::calcU(double tu, double tv) const
{
    return -tu * tv * tv + settings.F * (1.0 - tu);
}

double GrayScott::calcV(double tu, double tv) const
{
    return tu * tv * tv - (settings.F + settings.k) * tv;
}

double GrayScott::laplacian(int x, int y, int z, const std::vector<double> &s) const
{
    double ts = 0.0;
    ts += s[l2i(x - 1, y, z)];
    ts += s[l2i(x + 1, y, z)];
    ts += s[l2i(x, y - 1, z)];
    ts += s[l2i(x, y + 1, z)];
    ts += s[l2i(x, y, z - 1)];
    ts += s[l2i(x, y, z + 1)];
    ts += -6.0 * s[l2i(x, y, z)];

    return ts / 6.0;
}

void GrayScott::calc(const std::vector<double> &u, const std::vector<double> &v,
                     std::vector<double> &u2, std::vector<double> &v2)
{
    for (int z = 1, sizeZ = static_cast<int>(size_z); z < sizeZ + 1; z++)
    {
        for (int y = 1, sizeY = static_cast<int>(size_y); y < sizeY + 1; y++)
        {
            for (int x = 1, sizeX = static_cast<int>(size_x); x < sizeX + 1; x++)
            {
                const int i = l2i(x, y, z);
                double du = 0.0;
                double dv = 0.0;
                du = settings.Du * laplacian(x, y, z, u);
                dv = settings.Dv * laplacian(x, y, z, v);
                du += calcU(u[i], v[i]);
                dv += calcV(u[i], v[i]);
                du += settings.noise * uniform_dist(mt_gen);
                u2[i] = u[i] + du * settings.dt;
                v2[i] = v[i] + dv * settings.dt;
            }
        }
    }
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
    MPI_Type_vector(static_cast<int>(size_y) + 2, static_cast<int>(size_x),
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

void GrayScott::exchange_xy(std::vector<double> &local_data) const
{
    MPI_Status st;

    // Send XY face z=size_z to north and receive z=0 from south
    MPI_Sendrecv(&local_data[l2i(1, 0, static_cast<int>(size_z))], 1, xy_face_type, north, 1,
                 &local_data[l2i(1, 0, 0)], 1, xy_face_type, south, 1, cart_comm, &st);
    // Send XY face z=1 to south and receive z=size_z+1 from north
    MPI_Sendrecv(&local_data[l2i(1, 0, 1)], 1, xy_face_type, south, 1,
                 &local_data[l2i(1, 0, static_cast<int>(size_z + 1))], 1, xy_face_type, north, 1,
                 cart_comm, &st);
}

void GrayScott::exchange_xz(std::vector<double> &local_data) const
{
    MPI_Status st;

    // Send XZ face y=size_y to up and receive y=0 from down
    MPI_Sendrecv(&local_data[l2i(1, static_cast<int>(size_y), 1)], 1, xz_face_type, up, 2,
                 &local_data[l2i(1, 0, 1)], 1, xz_face_type, down, 2, cart_comm, &st);
    // Send XZ face y=1 to down and receive y=size_y+1 from up
    MPI_Sendrecv(&local_data[l2i(1, 1, 1)], 1, xz_face_type, down, 2,
                 &local_data[l2i(1, static_cast<int>(size_y + 1), 1)], 1, xz_face_type, up, 2,
                 cart_comm, &st);
}

void GrayScott::exchange_yz(std::vector<double> &local_data) const
{
    MPI_Status st;

    // Send YZ face x=size_x to east and receive x=0 from west
    MPI_Sendrecv(&local_data[l2i(static_cast<int>(size_x), 0, 0)], 1, yz_face_type, east, 3,
                 &local_data[l2i(0, 0, 0)], 1, yz_face_type, west, 3, cart_comm, &st);
    // Send YZ face x=1 to west and receive x=size_x+1 from east
    MPI_Sendrecv(&local_data[l2i(1, 0, 0)], 1, yz_face_type, west, 3,
                 &local_data[l2i(static_cast<int>(size_x + 1), 0, 0)], 1, yz_face_type, east, 3,
                 cart_comm, &st);
}

void GrayScott::exchange(std::vector<double> &u, std::vector<double> &v) const
{
    exchange_xy(u);
    exchange_xz(u);
    exchange_yz(u);

    exchange_xy(v);
    exchange_xz(v);
    exchange_yz(v);
}

void GrayScott::data_no_ghost_common(const std::vector<double> &data, double *data_no_ghost) const
{
    for (int z = 1, sizeZ = static_cast<int>(size_z); z < sizeZ + 1; z++)
    {
        for (int y = 1, sizeY = static_cast<int>(size_y); y < sizeY + 1; y++)
        {
            for (int x = 1, sizeX = static_cast<int>(size_x); x < sizeX + 1; x++)
            {
                data_no_ghost[(x - 1) + (y - 1) * size_x + (z - 1) * size_x * size_y] =
                    data[l2i(x, y, z)];
            }
        }
    }
}
