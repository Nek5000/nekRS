""" 
The present file contains runtime backend for using CPU Threads, and optionally 
CUDA.jl and AMDGPU.jl
"""
module Simulation

export init_domain, init_fields

import MPI
import Distributions

# from parent module
import ..Settings, ..MPICartDomain, ..Fields

# include functions for NVIDIA GPUs using CUDA.jl
include("Simulation_CUDA.jl")
# include functions for AMD GPUs using AMDGPU.jl
include("Simulation_AMDGPU.jl")

function init_domain(settings::Settings, comm::MPI.Comm)::MPICartDomain
    mcd = MPICartDomain()

    # set dims and Cartesian communicator
    mcd.dims = MPI.Dims_create(MPI.Comm_size(comm), mcd.dims)
    mcd.cart_comm = MPI.Cart_create(comm, mcd.dims)

    # set proc local coordinates in Cartesian communicator
    rank = MPI.Comm_rank(comm)
    mcd.coords = MPI.Cart_coords(mcd.cart_comm, rank)

    # set proc local mesh sizes
    mcd.proc_sizes = settings.L ./ mcd.dims

    for (i, coord) in enumerate(mcd.coords)
        if coord < settings.L % mcd.dims[i]
            mcd.proc_sizes[i] += 1
        end
    end

    # set proc local offsets
    for i in 1:3
        mcd.proc_offsets[i] = settings.L / mcd.dims[i] *
                              mcd.coords[i]
        +min(settings.L % mcd.dims[i], mcd.coords[i])
    end

    # get neighbors ranks
    mcd.proc_neighbors["west"],
    mcd.proc_neighbors["east"] = MPI.Cart_shift(mcd.cart_comm, 0, 1)
    mcd.proc_neighbors["down"],
    mcd.proc_neighbors["up"] = MPI.Cart_shift(mcd.cart_comm, 1, 1)
    mcd.proc_neighbors["south"],
    mcd.proc_neighbors["north"] = MPI.Cart_shift(mcd.cart_comm, 2, 1)

    return mcd
end

"""
Create and Initialize fields for either CPU, CUDA.jl, AMDGPU.jl backends
Multiple dispatch would direct to the appropriate overleaded function
"""
function init_fields(settings::Settings,
                     mcd::MPICartDomain, T)::Fields{T}
    lowercase_backend = lowercase(settings.backend)
    if lowercase_backend == "cuda"
        return _init_fields_cuda(settings, mcd, T)
    elseif lowercase_backend == "amdgpu"
        return _init_fields_amdgpu(settings, mcd, T)
    end
    # everything else would trigger the CPU threads backend
    return _init_fields_cpu(settings, mcd, T)
end

function _init_fields_cpu(settings::Settings,
                          mcd::MPICartDomain, T)::Fields{T}
    size_x = mcd.proc_sizes[1]
    size_y = mcd.proc_sizes[2]
    size_z = mcd.proc_sizes[3]

    # should be ones
    u = ones(T, size_x + 2, size_y + 2, size_z + 2)
    v = zeros(T, size_x + 2, size_y + 2, size_z + 2)

    u_temp = zeros(T, size_x + 2, size_y + 2, size_z + 2)
    v_temp = zeros(T, size_x + 2, size_y + 2, size_z + 2)

    function is_inside(x, y, z, offsets, sizes)::Bool
        if x < offsets[1] || x >= offsets[1] + sizes[1]
            return false
        end
        if y < offsets[2] || y >= offsets[2] + sizes[2]
            return false
        end
        if z < offsets[3] || z >= offsets[3] + sizes[3]
            return false
        end

        return true
    end

    d::Int64 = 6

    # global locations
    minL = Int64(settings.L / 2 - d)
    maxL = Int64(settings.L / 2 + d)

    xoff = mcd.proc_offsets[1]
    yoff = mcd.proc_offsets[2]
    zoff = mcd.proc_offsets[3]

    Threads.@threads for z in minL:maxL
        for y in minL:maxL
            for x in minL:maxL
                if !is_inside(x, y, z, mcd.proc_offsets, mcd.proc_sizes)
                    continue
                end

                # Julia is 1-index, like Fortran :)
                u[x - xoff + 2, y - yoff + 2, z - zoff + 2] = 0.25
                v[x - xoff + 2, y - yoff + 2, z - zoff + 2] = 0.33
            end
        end
    end

    xy_face_t, xz_face_t, yz_face_t = _get_mpi_faces(size_x, size_y, size_z, T)

    fields = Fields(u, v, u_temp, v_temp, xy_face_t, xz_face_t, yz_face_t)
    return fields
end

function iterate!(fields::Fields{T, N, Array{T, N}}, settings::Settings,
                  mcd::MPICartDomain) where {T, N}
    _exchange!(fields, mcd)
    # this function is the bottleneck
    _calculate!(fields, settings, mcd)

    # swap the names
    fields.u, fields.u_temp = fields.u_temp, fields.u
    fields.v, fields.v_temp = fields.v_temp, fields.v
end

function _get_mpi_faces(size_x, size_y, size_z, T)

    ## create a new type taking: count, block length, stride
    ## to interoperate with MPI for ghost cell exchange
    xy_face_t = MPI.Types.create_vector(size_y + 2, size_x, size_x + 2,
                                        MPI.Datatype(T))
    xz_face_t = MPI.Types.create_vector(size_z, size_x,
                                        (size_x + 2) * (size_y + 2),
                                        MPI.Datatype(T))
    yz_face_t = MPI.Types.create_vector((size_y + 2) * (size_z + 2), 1,
                                        size_x + 2, MPI.Datatype(T))
    MPI.Types.commit!(xy_face_t)
    MPI.Types.commit!(xz_face_t)
    MPI.Types.commit!(yz_face_t)

    return xy_face_t, xz_face_t, yz_face_t
end

function _exchange!(fields, mcd)
    """
    Send XY face z=size_z+1 to north and receive z=1 from south
    """
    function _exchange_xy!(var, size_z, data_type, rank1, rank2, comm)
        # to north
        send_buf = MPI.Buffer(@view(var[2, 1, size_z + 1]), 1, data_type)
        recv_buf = MPI.Buffer(@view(var[2, 1, 1]), 1, data_type)
        MPI.Sendrecv!(send_buf, recv_buf, comm, dest = rank1, source = rank2)

        # to south
        send_buf = MPI.Buffer(@view(var[2, 1, 2]), 1, data_type)
        recv_buf = MPI.Buffer(@view(var[2, 1, size_z + 2]), 1, data_type)
        MPI.Sendrecv!(send_buf, recv_buf, comm, dest = rank2, source = rank1)
    end

    """
    Send XZ face y=size_y+1 to up and receive y=1 from down
    """
    function _exchange_xz!(var, size_y, data_type, rank1, rank2, comm)
        # to up
        send_buf = MPI.Buffer(@view(var[2, size_y + 1, 2]), 1, data_type)
        recv_buf = MPI.Buffer(@view(var[2, 1, 2]), 1, data_type)
        MPI.Sendrecv!(send_buf, recv_buf, comm, dest = rank1, source = rank2)

        # to down
        send_buf = MPI.Buffer(@view(var[2, 2, 2]), 1, data_type)
        recv_buf = MPI.Buffer(@view(var[2, size_y + 2, 2]), 1, data_type)
        MPI.Sendrecv!(send_buf, recv_buf, comm, dest = rank2, source = rank1)
    end

    """
    Send YZ face x=size_x+2 to east and receive x=2 from west
    """
    function _exchange_yz!(var, size_x, data_type, rank1, rank2, comm)
        # to east
        send_buf = MPI.Buffer(@view(var[size_x + 1, 1, 1]), 1, data_type)
        recv_buf = MPI.Buffer(@view(var[1, 1, 1]), 1, data_type)
        MPI.Sendrecv!(send_buf, recv_buf, comm, dest = rank1, source = rank2)

        # to west
        send_buf = MPI.Buffer(@view(var[2, 1, 1]), 1, data_type)
        recv_buf = MPI.Buffer(@view(var[size_x + 2, 1, 1]), 1, data_type)
        MPI.Sendrecv!(send_buf, recv_buf, comm, dest = rank2, source = rank1)
    end

    # if already a CPU array, no need to copy, 
    # otherwise (device) copy to host. 
    u = typeof(fields.u) <: Array ? fields.u : Array(fields.u)
    v = typeof(fields.v) <: Array ? fields.v : Array(fields.v)

    for var in [u, v]
        _exchange_xy!(var, mcd.proc_sizes[3], fields.xy_face_t,
                      mcd.proc_neighbors["north"], mcd.proc_neighbors["south"],
                      mcd.cart_comm)

        _exchange_xz!(var, mcd.proc_sizes[2], fields.xz_face_t,
                      mcd.proc_neighbors["up"], mcd.proc_neighbors["down"],
                      mcd.cart_comm)

        _exchange_yz!(var, mcd.proc_sizes[1], fields.yz_face_t,
                      mcd.proc_neighbors["east"], mcd.proc_neighbors["west"],
                      mcd.cart_comm)
    end
end

function _calculate!(fields::Fields{T, N, Array{T, N}}, settings::Settings,
                     mcd::MPICartDomain) where {T, N}
    Du = convert(T, settings.Du)
    Dv = convert(T, settings.Dv)
    F = convert(T, settings.F)
    K = convert(T, settings.k)
    noise = convert(T, settings.noise)
    dt = convert(T, settings.dt)

    # loop through non-ghost cells, bounds are inclusive
    # @TODO: load balancing? option: a big linear loop
    # use @inbounds at the right for-loop level, avoid putting it at the top level
    Threads.@threads for k in 2:(mcd.proc_sizes[3] + 1)
        for j in 2:(mcd.proc_sizes[2] + 1)
            @inbounds for i in 2:(mcd.proc_sizes[1] + 1)
                u = fields.u[i, j, k]
                v = fields.v[i, j, k]

                # introduce a random disturbance on du
                du = Du * _laplacian(i, j, k, fields.u) - u * v^2 +
                     F * (1.0 - u) +
                     noise * rand(Distributions.Uniform(-1, 1))

                dv = Dv * _laplacian(i, j, k, fields.v) + u * v^2 -
                     (F + K) * v

                # advance the next step
                fields.u_temp[i, j, k] = u + du * dt
                fields.v_temp[i, j, k] = v + dv * dt
            end
        end
    end
end

"""
   7-point stencil around the cell, 
   this is equally a host and a device function!
"""
function _laplacian(i, j, k, var)
    @inbounds l = var[i - 1, j, k] + var[i + 1, j, k] + var[i, j - 1, k] +
                  var[i, j + 1, k] + var[i, j, k - 1] + var[i, j, k + 1] -
                  6.0 * var[i, j, k]
    return l / 6.0
end

function get_fields(fields::Fields{T, N, Array{T, N}}) where {T, N}
    @inbounds begin
        u_no_ghost = fields.u[(begin + 1):(end - 1), (begin + 1):(end - 1),
                              (begin + 1):(end - 1)]
        v_no_ghost = fields.v[(begin + 1):(end - 1), (begin + 1):(end - 1),
                              (begin + 1):(end - 1)]
    end
    return u_no_ghost, v_no_ghost
end

end # module 
