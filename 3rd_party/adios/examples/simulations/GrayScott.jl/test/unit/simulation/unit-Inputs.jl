
import Test: @testset, @test, @test_throws
import GrayScott: Inputs

# Unfortunately due to MPI being a Singleton, single MPI.Init()
# these unit tests don't run as independent files

@testset "unit-Inputs.get_settings" begin
    config_file = joinpath(dirname(Base.active_project()), "examples",
                           "settings-files.json")
    Inputs.get_settings([config_file], MPI.COMM_WORLD)

    @test_throws(ArgumentError,
                 Inputs.get_settings(["hello.nojson"], MPI.COMM_WORLD))
end
