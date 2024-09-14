
import Test: @testset, @test, @test_throws

import ADIOS2

# import submodule
import GrayScott: IO
# import types
import GrayScott: Settings, MPICartDomain, Fields

@testset "unit-IO.init" begin
    settings = Settings()
    mpi_cart_domain = MPICartDomain()
    fields = Simulation.init_fields(settings, mpi_cart_domain, Float32)

    @test eltype(fields.u) == Float32
    @test eltype(fields.v) == Float32

    stream = IO.init(settings, mpi_cart_domain, fields)

    @test ADIOS2.name(stream.engine) == "foo.bp"
    IO.close!(stream)

    # @TODO: needs to be done from rank==0 only 
    # Base.Filesystem.rm("foo.bp", force = true, recursive = true)
end

@testset "unit-IO.write" begin
    settings = Settings()
    settings.L = 6
    mpi_cart_domain = Simulation.init_domain(settings, MPI.COMM_WORLD)
    fields = Simulation.init_fields(settings, mpi_cart_domain, Float32)
    stream = IO.init(settings, mpi_cart_domain, fields)
    IO.write_step!(stream, Int32(0), fields)
    IO.close!(stream)
end
