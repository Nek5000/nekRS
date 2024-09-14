"""
GrayScott.jl is a Simulation and Analysis parallel framework for solving the 
Gray-Scott 3D diffusion reaction system of equations of two variables U and V on 
a regular Cartesian mesh.

The bp output files can be visualized with ParaView.
"""
module GrayScott

import MPI, ADIOS2

# contains relevant data containers "structs" for Input, Domain and Fields
include(joinpath("simulation", "Structs.jl"))

# contains helper functions for general use
include(joinpath("helper", "Helper.jl"))
import .Helper

# initializes inputs from configuration file
include(joinpath("simulation", "Inputs.jl"))
import .Inputs

# manages the simulation computation
include(joinpath("simulation", "Simulation.jl"))
import .Simulation

# manages the I/O
include(joinpath("simulation", "IO.jl"))
import .IO

function julia_main()::Cint
    try
        main(ARGS)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function main(args::Vector{String})::Int32
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    # a data struct that holds settings data from config_file in args
    # example config file: ../examples/settings-files.json
    settings = Inputs.get_settings(args, comm)

    # initialize MPI Cartesian Domain and Communicator
    mpi_cart_domain = Simulation.init_domain(settings, comm)

    # initialize fields
    fields = Simulation.init_fields(settings,
                                    mpi_cart_domain,
                                    Helper.get_type(settings.precision))

    # initialize IOStream struct holding ADIOS-2 components for parallel I/O
    stream = IO.init(settings, mpi_cart_domain, fields)

    restart_step::Int32 = 0
    # @TODO: checkpoint-restart 
    step::Int32 = restart_step

    while step < settings.steps
        Simulation.iterate!(fields, settings, mpi_cart_domain)
        step += 1

        if step % settings.plotgap == 0
            if rank == 0
                println("Simulation at step ", step, " writing output step ",
                        step / settings.plotgap)
            end

            IO.write_step!(stream, step, fields)
        end
    end

    IO.close!(stream)

    # Debugging session or Julia REPL session, not needed overall as it would be 
    # called when the program ends
    if !isinteractive()
        MPI.Finalize()
    end

    return 0
end

end # module GrayScott
