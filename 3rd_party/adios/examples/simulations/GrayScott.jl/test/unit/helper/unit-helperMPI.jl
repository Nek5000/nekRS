
import Test: @testset, @test, @test_throws
import GrayScott: Helper

@testset "unit-Helper.bcast_file_contents" begin
    config_file = joinpath(dirname(Base.active_project()), "examples",
                           "settings-files.json")

    file_contents = Helper.bcast_file_contents(config_file, MPI.COMM_WORLD)
    file_contents_expected = String(read(open(config_file, "r")))
    @test file_contents == file_contents_expected
end