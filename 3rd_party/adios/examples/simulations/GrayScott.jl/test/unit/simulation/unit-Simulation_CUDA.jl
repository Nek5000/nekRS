
import Test: @testset, @test, @test_throws
# import submodule
import GrayScott: Simulation
# import types
import GrayScott: Settings, MPICartDomain, Fields

@testset "unit-Simulation.init_fields-cuda" begin
    function test_init_cuda(L)
        settings = Settings()
        settings.L = L
        mpi_cart_domain = Simulation.init_domain(settings, MPI.COMM_WORLD)

        fields = Simulation.init_fields(settings, mpi_cart_domain, Float32)

        settings.backend = "CUDA"
        fields_cuda = Simulation.init_fields(settings, mpi_cart_domain, Float32)

        @test fields.u ≈ Array(fields_cuda.u)
        @test fields.v ≈ Array(fields_cuda.v)
    end

    test_init_cuda(8)
    test_init_cuda(16)
    test_init_cuda(32)
    test_init_cuda(64)
    test_init_cuda(128)
end