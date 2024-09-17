
"""
Settings carry the settings from the simulation config file (json or yaml formats)

Using Base.@kwdef macro for easy defaults and enable keyword arguments
Settings(Du = 0.2, noise = 0.2) 
See:
https://discourse.julialang.org/t/default-value-of-some-fields-in-a-mutable-struct/33408/24?u=williamfgc
"""
Base.@kwdef mutable struct Settings
    L::Int64 = 128
    steps::Int32 = 20000
    plotgap::Int32 = 200
    F::Float64 = 0.04
    k::Float64 = 0
    dt::Float64 = 0.2
    Du::Float64 = 0.05
    Dv::Float64 = 0.1
    noise::Float64 = 0.0
    output::String = "foo.bp"
    checkpoint::Bool = false
    checkpoint_freq::Int32 = 2000
    checkpoint_output::String = "ckpt.bp"
    restart::Bool = false
    restart_input::String = "ckpt.bp"
    adios_config::String = "adios2.yaml"
    adios_span::Bool = false
    adios_memory_selection::Bool = false
    mesh_type::String = "image"
    precision::String = "Float64"
    backend::String = "CPU"
end

SettingsKeys = Set{String}([
                               "L",
                               "steps",
                               "plotgap",
                               "F",
                               "k",
                               "dt",
                               "Du",
                               "Dv",
                               "noise",
                               "output",
                               "checkpoint",
                               "checkpoint_freq",
                               "checkpoint_output",
                               "restart",
                               "restart_input",
                               "adios_config",
                               "adios_span",
                               "adios_memory_selection",
                               "mesh_type",
                               "precision",
                               "backend",
                           ])

Base.@kwdef mutable struct MPICartDomain
    cart_comm::MPI.Comm = MPI.COMM_NULL

    # Cartesian communicator info
    # Could used StaticArrays.jl?
    # start dims with zeros
    dims::Vector{Int32} = zeros(Int32, 3)
    coords::Vector{Int32} = zeros(Int32, 3)

    # local process mesh sizes and offsets in Cartesian domain info, using defaults
    proc_sizes::Vector{Int64} = [128, 128, 128]
    proc_offsets::Vector{Int64} = [1, 1, 1]

    # couldn't use NamedTuples as struct is mutable
    proc_neighbors = Dict{String, Int32}("west" => -1, "east" => -1, "up" => -1,
                                         "down" => -1, "north" => -1,
                                         "south" => -1)
end

"""
Carry the physical field outputs: u and v
Using AbstractArray to allow for Array, CuArray and ROCArray
"""
mutable struct Fields{T, N, A <: AbstractArray{T, N}}
    u::A
    v::A
    u_temp::A
    v_temp::A
    # MPI Datatypes for halo exchange MPI.Datatype(T)
    xy_face_t::MPI.Datatype
    xz_face_t::MPI.Datatype
    yz_face_t::MPI.Datatype
end

"""
Carry the I/O information for outputs
"""
struct IOStream
    adios::ADIOS2.Adios
    io::ADIOS2.AIO
    engine::ADIOS2.Engine
    var_step::ADIOS2.Variable
    var_U::ADIOS2.Variable
    var_V::ADIOS2.Variable
end
