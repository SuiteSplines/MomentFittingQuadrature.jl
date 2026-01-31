using Test

tests = [
    "algoim",
    "bernstein",
    "weights",
    "stencils"
]

@testset "MomentFittingQuadrature" begin
    for t in tests
        fp = joinpath(dirname(@__FILE__), "$t.jl")
        println("$fp ...")
        include(fp)
    end
end # @testset
