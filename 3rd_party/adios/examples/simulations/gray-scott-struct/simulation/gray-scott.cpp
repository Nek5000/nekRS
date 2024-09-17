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
#include <stddef.h>  // offsetof
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
    exchange(d);
    calc(d, d2);
    d.swap(d2);
}

void GrayScott::restart(std::vector<MemLayout> &d_in)
{
    auto expected_len = (size_x + 2) * (size_y + 2) * (size_z + 2);
    if (d_in.size() == expected_len)
    {
        d = d_in;
    }
    else
    {
        throw std::runtime_error("Restart with incompatible array size, expected " +
                                 std::to_string(expected_len) + " got " +
                                 std::to_string(d_in.size()) + " elements");
    }
}

const std::vector<GrayScott::MemLayout> &GrayScott::d_ghost() const { return d; }

std::vector<GrayScott::MemLayout> GrayScott::d_noghost() const { return data_noghost(d); }

void GrayScott::d_noghost(GrayScott::MemLayout *u_no_ghost) const { data_noghost(d, u_no_ghost); }

std::vector<GrayScott::MemLayout>
GrayScott::data_noghost(const std::vector<GrayScott::MemLayout> &data) const
{
    std::vector<GrayScott::MemLayout> buf(size_x * size_y * size_z);
    data_no_ghost_common(data, buf.data());
    return buf;
}

void GrayScott::data_noghost(const std::vector<GrayScott::MemLayout> &data,
                             GrayScott::MemLayout *data_no_ghost) const
{
    data_no_ghost_common(data, data_no_ghost);
}

void GrayScott::init_field()
{
    const int V = static_cast<int>((size_x + 2) * (size_y + 2) * (size_z + 2));
    MemLayout initvalue = {1.0, 0.0};
    d.resize(V, initvalue);
    d2.resize(V, initvalue);

    const int dd = 6;
    const auto settingsL = static_cast<int>(settings.L);
    for (int z = settingsL / 2 - dd; z < settingsL / 2 + dd; z++)
    {
        for (int y = settingsL / 2 - dd; y < settingsL / 2 + dd; y++)
        {
            for (int x = settingsL / 2 - dd; x < settingsL / 2 + dd; x++)
            {
                if (!is_inside(x, y, z))
                    continue;
                int i = g2i(x, y, z);
                d[i].u = 0.25;
                d[i].v = 0.33;
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

GrayScott::MemLayout GrayScott::laplacian(int x, int y, int z,
                                          const std::vector<MemLayout> &s) const
{
    MemLayout dd = {0.0, 0.0};
    {
        double &ts = dd.u;
        ts += s[l2i(x - 1, y, z)].u;
        ts += s[l2i(x + 1, y, z)].u;
        ts += s[l2i(x, y - 1, z)].u;
        ts += s[l2i(x, y + 1, z)].u;
        ts += s[l2i(x, y, z - 1)].u;
        ts += s[l2i(x, y, z + 1)].u;
        ts += -6.0 * s[l2i(x, y, z)].u;
        ts = ts / 6.0;
    }

    {
        double &ts = dd.v;
        ts += s[l2i(x - 1, y, z)].v;
        ts += s[l2i(x + 1, y, z)].v;
        ts += s[l2i(x, y - 1, z)].v;
        ts += s[l2i(x, y + 1, z)].v;
        ts += s[l2i(x, y, z - 1)].v;
        ts += s[l2i(x, y, z + 1)].v;
        ts += -6.0 * s[l2i(x, y, z)].v;
        ts = ts / 6.0;
    }
    return dd;
}

void GrayScott::calc(const std::vector<MemLayout> &d, std::vector<MemLayout> &d2)
{
    for (int z = 1, sizeZ = static_cast<int>(size_z); z < sizeZ + 1; z++)
    {
        for (int y = 1, sizeY = static_cast<int>(size_y); y < sizeY + 1; y++)
        {
            for (int x = 1, sizeX = static_cast<int>(size_x); x < sizeX + 1; x++)
            {
                const int i = l2i(x, y, z);
                MemLayout dd = laplacian(x, y, z, d);
                dd.u = settings.Du * dd.u;
                dd.v = settings.Dv * dd.v;
                dd.u += calcU(d[i].u, d[i].v);
                dd.v += calcV(d[i].u, d[i].v);
                dd.u += settings.noise * uniform_dist(mt_gen);
                d2[i].u = d[i].u + dd.u * settings.dt;
                d2[i].v = d[i].v + dd.v * settings.dt;
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

    // MPI type for Memory Layout
    int count = 2;
    int array_of_blocklengths[] = {1, 1};
    MPI_Aint array_of_displacements[] = {offsetof(MemLayout, u), offsetof(MemLayout, v)};
    MPI_Datatype array_of_types[] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types,
                           &memlayout_type);
    MPI_Type_commit(&memlayout_type);

    // XY faces: size_x * (size_y + 2)
    MPI_Type_vector(static_cast<int>(size_y + 2), static_cast<int>(size_x),
                    static_cast<int>(size_x + 2), memlayout_type, &xy_face_type);
    MPI_Type_commit(&xy_face_type);

    // XZ faces: size_x * size_z
    MPI_Type_vector(static_cast<int>(size_z), static_cast<int>(size_x),
                    static_cast<int>((size_x + 2) * (size_y + 2)), memlayout_type, &xz_face_type);
    MPI_Type_commit(&xz_face_type);

    // YZ faces: (size_y + 2) * (size_z + 2)
    MPI_Type_vector(static_cast<int>((size_y + 2) * (size_z + 2)), 1, static_cast<int>(size_x + 2),
                    memlayout_type, &yz_face_type);
    MPI_Type_commit(&yz_face_type);
}

void GrayScott::exchange_xy(std::vector<MemLayout> &local_data) const
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

void GrayScott::exchange_xz(std::vector<MemLayout> &local_data) const
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

void GrayScott::exchange_yz(std::vector<MemLayout> &local_data) const
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

void GrayScott::exchange(std::vector<MemLayout> &d) const
{
    exchange_xy(d);
    exchange_xz(d);
    exchange_yz(d);
}

void GrayScott::data_no_ghost_common(const std::vector<MemLayout> &data,
                                     MemLayout *data_no_ghost) const
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
