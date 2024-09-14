
import MPI
import ArgParse
import ADIOS2

function _epsilon(d::T)::Bool where {T <: Number}
    return (d < 1.0e-20)
end

"""
Return 2 arrays pdf and bins of a 2D slice
"""
function _compute_pdf(data::Array{T, 3}, shape, count, nbins, min::T,
                      max::T) where {T}
    pdf = Array{T, 2}(undef, nbins, count)
    bins = Array{T, 1}(undef, nbins)
    bin_width = (max - min) / nbins

    for i in 1:nbins
        bins[i] = min + (i - 1) * bin_width
    end

    slice_size = shape[2] * shape[3]
    # special case: only one bin or small window
    if nbins == 1 || _epsilon(max - min) || _epsilon(bin_width)
        fill!(pdf, slice_size)
        return pdf, bins
    end

    # Calculate a PDF for 'nbins' bins for values between 'min' and 'max'

    for c in 1:count
        for j in 1:shape[2]
            for k in 1:shape[3]
                value = data[k, j, c]

                if value > max || value < min
                    @show value, " is outside [", min, max "]"
                end
                bin = floor((value - min) / bin_width)

                if bin == nbins
                    bin = nbins - 1
                end
                pdf[bin, c] += 1
            end
        end
    end
end

function _parse_arguments(args)
    s = ArgParse.ArgParseSettings(description = "gray-scott workflow pdf generator, Julia version")

    #  @add_arg_table! s begin
    #       "--opt1"               # an option (will take an argument)
    #       "--opt2", "-o"         # another option, with short form
    #       "arg1"                 # a positional argument
    #   end

    ArgParse.@add_arg_table! s begin
        "input"
        help = "Name of the input file handle for reading data"
        arg_type = String
        required = true
        "output"
        help = "Name of the output file to which data must be written"
        arg_type = String
        required = true
        "N"
        help = "Number of bins for the PDF calculation, default = 1000"
        arg_type = Int64
        required = false
        default = 1000
        "output_inputdata"
        help = "YES will write the original variables besides the analysis results"
        arg_type = Bool
        required = false
        default = false
    end

    # parse_args return a dictionary with key/value for arguments
    parsed_arguments = ArgParse.parse_args(args, s)
    return parsed_arguments
end

function _read_data_write_pdf(inputs, comm)
    in_filename = inputs["input"]
    out_filename = inputs["output"]
    nbins = inputs["N"]
    write_inputvars = inputs["output_inputdata"]

    adios = ADIOS2.adios_init_mpi("adios2.xml", comm)

    reader_io = ADIOS2.declare_io(adios, "SimulationOutput")
    writer_io = ADIOS2.declare_io(adios, "PDFAnalysisOutput")

    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    if rank == 0
        println("PDF analysis reads from Simulation using engine type:  ",
                ADIOS2.engine_type(reader_io))
        println("PDF analysis writes using engine type:  ",
                ADIOS2.engine_type(writer_io))
    end

    # Engines for reading and writing
    reader = ADIOS2.open(reader_io, inputs["input"], ADIOS2.mode_read)
    writer = ADIOS2.open(writer_io, inputs["output"], ADIOS2.mode_write)

    # break inside if needed
    while true

        # timeout is 10 seconds
        read_status = ADIOS2.begin_step(reader, ADIOS2.step_mode_read, 10)

        if read_status == ADIOS2.step_status_not_ready
            # sleep in seconds, minimum is one milisecond = 0.001
            sleep(1)
            continue
        else if read_status != ADIOS2.step_status_ok
            break
        end

        step_sim_out = ADIOS2.current_step(reader)
        var_U = ADIOS2.inquire_variable(reader_io, "U")
        var_V = ADIOS2.inquire_variable(reader_io, "V")
        var_step = ADIOS2.inquire_variable(reader_io, "step")

        shape = ADIOS2.shape(var_U)

        # Split in the slowest dimension
        count_z = shape[3] / size
        start_z = count_z * rank

        # Last process needs to read all the rest
        if rank == size-1
            count_z = shape[3] - count_z * (size-1)
        end

        # missing set_selection
        start = ( 0,0,start_z)
        count = (shape[1], shape[2], count_z)
        ADIOS2.set_selection(var_U, start, count)
        ADIOS2.set_selection(var_V, start, count)

        # Calculate

    end
end

function main(args)
    MPI.Init()

    inputs = _parse_arguments(args)
    comm = MPI.COMM_WORLD
    _read_data_write_pdf(inputs, comm)

    MPI.Finalize()
end
