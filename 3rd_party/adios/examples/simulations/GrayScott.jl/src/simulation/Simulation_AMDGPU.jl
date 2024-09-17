
import AMDGPU

function _init_fields_amdgpu(settings::Settings, mcd::MPICartDomain,
                             T)::Fields{T, 3, <:AMDGPU.ROCArray{T, 3}}
    size_x = mcd.proc_sizes[1]
    size_y = mcd.proc_sizes[2]
    size_z = mcd.proc_sizes[3]

    # should be ones
    u = AMDGPU.ones(T, size_x + 2, size_y + 2, size_z + 2)
    v = AMDGPU.zeros(T, size_x + 2, size_y + 2, size_z + 2)

    u_temp = AMDGPU.zeros(T, size_x + 2, size_y + 2, size_z + 2)
    v_temp = AMDGPU.zeros(T, size_x + 2, size_y + 2, size_z + 2)

    roc_offsets = AMDGPU.ROCArray(mcd.proc_offsets)
    roc_sizes = AMDGPU.ROCArray(mcd.proc_sizes)

    d::Int64 = 6
    minL = Int64(settings.L / 2 - d)
    maxL = Int64(settings.L / 2 + d)

    # @TODO: get ideal grid size and threads
    threads = (16, 16)
    # grid size must be the total number of threads of each direction
    grid = (settings.L, settings.L)

    AMDGPU.wait(AMDGPU.@roc groupsize=threads gridsize=grid _populate_amdgpu!(u,
                                                                              v,
                                                                              roc_offsets,
                                                                              roc_sizes,
                                                                              minL,
                                                                              maxL))

    xy_face_t, xz_face_t, yz_face_t = _get_mpi_faces(size_x, size_y, size_z, T)

    fields = Fields(u, v, u_temp, v_temp, xy_face_t, xz_face_t, yz_face_t)
    return fields
end

function iterate!(fields::Fields{T, N, <:AMDGPU.ROCArray{T, N}},
                  settings::Settings,
                  mcd::MPICartDomain) where {T, N}
    _exchange!(fields, mcd)
    # this function is the bottleneck
    _calculate!(fields, settings, mcd)

    # swap the names
    fields.u, fields.u_temp = fields.u_temp, fields.u
    fields.v, fields.v_temp = fields.v_temp, fields.v
end

function _populate_amdgpu!(u, v, offsets, sizes, minL, maxL)
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

    # local coordinates (this are 1-index already)
    lz = (AMDGPU.workgroupIdx().x - Int32(1)) * AMDGPU.workgroupDim().x +
         AMDGPU.workitemIdx().x
    ly = (AMDGPU.workgroupIdx().y - Int32(1)) * AMDGPU.workgroupDim().y +
         AMDGPU.workitemIdx().y

    # This check might not be needed
    if lz <= size(u, 3) && ly <= size(u, 2)

        # get global coordinates
        z = lz + offsets[3] - 1
        y = ly + offsets[2] - 1

        if z >= minL && z <= maxL && y >= minL && y <= maxL
            xoff = offsets[1]

            for x in minL:maxL
                # check if global coordinates for initialization are inside the region
                if !is_inside(x, y, z, offsets, sizes)
                    continue
                end

                # Julia is 1-index, like Fortran :)
                u[x - xoff + 2, ly + 1, lz + 1] = 0.25
                v[x - xoff + 2, ly + 1, lz + 1] = 0.33
            end
        end
    end
end

function _calculate!(fields::Fields{T, N, <:AMDGPU.ROCArray{T, N}},
                     settings::Settings,
                     mcd::MPICartDomain) where {T, N}
    Du = convert(T, settings.Du)
    Dv = convert(T, settings.Dv)
    F = convert(T, settings.F)
    K = convert(T, settings.k)
    noise = convert(T, settings.noise)
    dt = convert(T, settings.dt)

    roc_sizes = AMDGPU.ROCArray(mcd.proc_sizes)

    threads = (16, 16)
    blocks = (settings.L, settings.L)

    AMDGPU.wait(AMDGPU.@roc groupsize=threads gridsize=grid _calculate_kernel_amdgpu!(fields.u,
                                                                                      fields.v,
                                                                                      fields.u_temp,
                                                                                      fields.v_temp,
                                                                                      roc_sizes,
                                                                                      Du,
                                                                                      Dv,
                                                                                      F,
                                                                                      K,
                                                                                      noise,
                                                                                      dt))
end

function _calculate_kernel_amdgpu!(u, v, u_temp, v_temp, sizes, Du, Dv, F, K,
                                   noise, dt)

    # local coordinates (this are 1-index already)
    k = (AMDGPU.workgroupIdx().x - Int32(1)) * AMDGPU.workgroupDim().x +
        AMDGPU.workitemIdx().x
    j = (AMDGPU.workgroupIdx().y - Int32(1)) * AMDGPU.workgroupDim().y +
        AMDGPU.workitemIdx().y

    # loop through non-ghost cells
    if k >= 2 && k <= sizes[3] + 1 && j >= 2 && j <= sizes[2] + 1
        # bounds are inclusive
        for i in 2:(sizes[1] + 1)
            u_ijk = u[i, j, k]
            v_ijk = v[i, j, k]

            du = Du * _laplacian(i, j, k, u) - u_ijk * v_ijk^2 +
                 F * (1.0 - u_ijk)
            # + noise * AMDGPU.rand(eltype(u))
            # WIP in AMDGPU.jl, works with CUDA.jl
            # + rand(Distributions.Uniform(-1, 1))

            dv = Dv * _laplacian(i, j, k, v) + u_ijk * v_ijk^2 -
                 (F + K) * v_ijk

            # advance the next step
            u_temp[i, j, k] = u_ijk + du * dt
            v_temp[i, j, k] = v_ijk + dv * dt
        end
    end
end

function get_fields(fields::Fields{T, N, <:AMDGPU.ROCArray{T, N}}) where {T, N}
    u = Array(fields.u)
    u_no_ghost = u[(begin + 1):(end - 1), (begin + 1):(end - 1),
                   (begin + 1):(end - 1)]

    v = Array(fields.v)
    v_no_ghost = v[(begin + 1):(end - 1), (begin + 1):(end - 1),
                   (begin + 1):(end - 1)]
    return u_no_ghost, v_no_ghost
end
