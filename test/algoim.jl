using MomentFittingQuadrature, IgaFormation, SortedSequences

@testset "2d element corners" begin
    e = Element(Interval(1.0,2.0) ⨱ Interval(3.0,4.0))
    @test MomentFittingQuadrature.llcorner(e) == [1.0,3.0]
    @test MomentFittingQuadrature.urcorner(e) == [2.0,4.0]
end

@testset "3d element corners" begin
    e = Element(Interval(1.0,2.0) ⨱ Interval(3.0,4.0) ⨱ Interval(5.0,6.0))
    @test MomentFittingQuadrature.llcorner(e) == [1.0,3.0,5.0]
    @test MomentFittingQuadrature.urcorner(e) == [2.0,4.0,6.0]
end

@testset "Algoim element quadrature rule" begin
    e = Element(Interval(-1.0,1.0) ⨱ Interval(-1.0,1.0))
    ϕ(x) = sqrt(x[1]^2 + x[2]^2) - 1.0
    ∇ϕ(x) = [x[1]; x[2]] / sqrt(x[1]^2 + x[2]^2)

    x, w = algoim_quad_data(ϕ, ∇ϕ, e; degree=16, phase=1)
    @test isapprox(sum(w), 4-π; atol=10e-5)

    x, w = algoim_quad_data(ϕ, ∇ϕ, e; degree=16, phase=-1)
    @test isapprox(sum(w), π; atol=10e-5)

    x, w = algoim_quad_data(ϕ, ∇ϕ, e; degree=16, phase=0)
    @test isapprox(sum(w), 2π; atol=10e-4)
end

@testset "Algoim element quadrature rule using ForwardDiff" begin
    e = Element(Interval(-1.0,1.0) ⨱ Interval(-1.0,1.0))
    ϕ(x) = sqrt(x[1]^2 + x[2]^2) - 1.0

    x, w = algoim_quad_data(ϕ, e; degree=16, phase=1)
    @test isapprox(sum(w), 4-π; atol=10e-5)

    x, w = algoim_quad_data(ϕ, e; degree=16, phase=-1)
    @test isapprox(sum(w), π; atol=10e-5)

    x, w = algoim_quad_data(ϕ, e; degree=16, phase=0)
    @test isapprox(sum(w), 2π; atol=10e-4)
end
