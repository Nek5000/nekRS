module IO

export init, write_step!

# import external module
import ADIOS2

# internal submodule
import ..Simulation
# types from parent module types in Structs.jl
import ..Settings, ..MPICartDomain, ..Fields, ..IOStream

function init(settings::Settings, mcd::MPICartDomain,
              fields::Fields{T}) where {T}

    # initialize adios MPI using the cartesian communicator
    adios = ADIOS2.adios_init_mpi(mcd.cart_comm)
    io = ADIOS2.declare_io(adios, "SimulationOutput")
    # @TODO: implement ADIOS2.set_engine in ADIOS2.jl
    engine = ADIOS2.open(io, settings.output, ADIOS2.mode_write)

    # store simulation run provenance as attributes
    ADIOS2.define_attribute(io, "F", settings.F)
    ADIOS2.define_attribute(io, "k", settings.k)
    ADIOS2.define_attribute(io, "dt", settings.dt)
    ADIOS2.define_attribute(io, "Du", settings.Du)
    ADIOS2.define_attribute(io, "Dv", settings.Dv)
    ADIOS2.define_attribute(io, "noise", settings.noise)

    _add_visualization_schemas(io, settings.L)

    # ADIOS2 requires tuples for the dimensions
    # define global variables u and v
    shape = (settings.L, settings.L, settings.L)
    start = Tuple(mcd.proc_offsets)
    count = Tuple(mcd.proc_sizes)

    var_step = ADIOS2.define_variable(io, "step", Int32)
    var_U = ADIOS2.define_variable(io, "U", T, shape, start, count)
    var_V = ADIOS2.define_variable(io, "V", T, shape, start, count)

    return IOStream(adios, io, engine, var_step, var_U, var_V)
end

function write_step!(stream::IOStream, step::Int32, fields::Fields{T}) where {T}

    # this creates temporaries similar to Fortran
    u_no_ghost, v_no_ghost = Simulation.get_fields(fields)

    # writer engine
    w = stream.engine

    ADIOS2.begin_step(w)
    ADIOS2.put!(w, stream.var_step, step)
    ADIOS2.put!(w, stream.var_U, u_no_ghost)
    ADIOS2.put!(w, stream.var_V, v_no_ghost)
    ADIOS2.end_step(w)
end

function close!(stream::IOStream)
    ADIOS2.close(stream.engine)
    ADIOS2.adios_finalize(stream.adios)
end

function _add_visualization_schemas(io, length)

    # Fides schema
    ADIOS2.define_attribute(io, "Fides_Data_Model", "uniform")
    ADIOS2.define_attribute_array(io, "Fides_Origin", [0.0, 0.0, 0.0])
    ADIOS2.define_attribute_array(io, "Fides_Spacing", [0.1, 0.1, 0.1])
    ADIOS2.define_attribute(io, "Fides_Dimension_Variable", "U")
    ADIOS2.define_attribute_array(io, "Fides_Variable_List", ["U", "V"])
    ADIOS2.define_attribute_array(io, "Fides_Variable_Associations",
                                  ["points", "points"])

    # VTX schema
    # string concatenation uses *, ^ is for repetition 
    # if length = 64
    # extent = "0 64 0 64 0 64 "
    extent = ("0 " * string(length) * " ")^3
    extent = rstrip(extent)

    # deactive code formatting using JuliaFormatter.jl 
    # raw strings: " must be escaped with \"
    #! format: off
    vtx_schema = raw"
        <?xml version=\"1.0\"?>
        <VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">
          <ImageData WholeExtent=\"" * extent * raw"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">
            <Piece Extent=\"" * extent * raw"\">
              <CellData Scalars=\"U\">
                <DataArray Name=\"U\" />
                <DataArray Name=\"V\" />
                <DataArray Name=\"TIME\">
                  step
                </DataArray>
              </CellData>
            </Piece>
          </ImageData>
        </VTKFile>"
    #! format: on
    # reactive code formatting using JuliaFormatter.jl
    ADIOS2.define_attribute(io, "vtk.xml", vtx_schema)
end

end
