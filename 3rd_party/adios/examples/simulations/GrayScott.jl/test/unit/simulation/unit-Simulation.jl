
import Test: @testset, @test, @test_throws
# import submodule
import GrayScott: Simulation
# import types
import GrayScott: Settings, MPICartDomain, Fields

# Unfortunately due to MPI being a Singleton, single MPI.Init()
# these unit tests don't run as independent files

@testset "unit-Simulation.init" begin
    settings = Settings()
    mpi_cart_domain = Simulation.init_domain(settings, MPI.COMM_WORLD)
    fields = Simulation.init_fields(settings, mpi_cart_domain, Float32)

    @test typeof(fields) == Fields{Float32, 3, Array{Float32, 3}}
end

@testset "unit-Simulation.iterate" begin
    settings = Settings()
    settings.L = 2
    mpi_cart_domain = Simulation.init_domain(settings, MPI.COMM_WORLD)
    fields = Simulation.init_fields(settings, mpi_cart_domain, Float32)

    Simulation.iterate!(fields, settings, mpi_cart_domain)

    if verbose
        sleep(0.01)
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        @show rank, fields.v
    end
end

#end
