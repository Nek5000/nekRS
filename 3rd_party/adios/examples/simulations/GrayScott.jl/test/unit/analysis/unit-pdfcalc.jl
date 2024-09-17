
import Test: @testset, @test, @test_throws

include(joinpath(dirname(Base.active_project()), "src", "analysis",
                 "pdfcalc.jl"))

@testset "unit-analysis.pdfcalc._parse_args" begin
    inputs = _parse_arguments(["foo.bp", "bar.bp", "1500"])
    @test inputs["input"] == "foo.bp"
    @test inputs["output"] == "bar.bp"
    @test inputs["N"] == 1500
    @test inputs["output_inputdata"] == false

    inputs = _parse_arguments(["input.bp", "output.bp", "1000", "true"])
    @test inputs["input"] == "input.bp"
    @test inputs["output"] == "output.bp"
    @test inputs["N"] == 1000
    @test inputs["output_inputdata"] == true
end
