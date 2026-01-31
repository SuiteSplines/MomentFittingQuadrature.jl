@testset "Reshape elementwise data to patch data" begin
    patch_data = [
        1   7  13  19  25  31
        2   8  14  20  26  32
        3   9  15  21  27  33
        4  10  16  22  28  34
        5  11  17  23  29  35
        6  12  18  24  30  36
    ] # 3×2 blocks in column major order map to element data

    element_data = [
        1   4  13  16  25  28
        2   5  14  17  26  29
        3   6  15  18  27  30
        7  10  19  22  31  34
        8  11  20  23  32  35
        9  12  21  24  33  36
    ] # columns map to 3×2 blocks in patch data

    npoints = (3,2)
    nelem = (2, 3)

    # test reshaping
    @test weights_patch_reshape(element_data, npoints, nelem) ≈ patch_data
end

@testset "Moment fitting reference data and cache" begin
    A, b = MomentFittingQuadrature.moment_fitting_cache(; npoints=(4,2))

    # test dimensions
    @test size(A,1) == size(A,2)
    @test size(A,1) == size(b,1)

    # test zero b
    @test all(b .== 0)
end

@testset "Moment fitting on single element" begin
    # element
    e = Element(Interval(0.0,2.0) ⨱ Interval(0.0,2.0))

    # levelset
    ϕ(x) = sqrt(x[1]^2 + x[2]^2) - 2.0

    # number of quadrature points
    npoints = (3,3)

    # number of element subdivisions for moment integration
    algoim_nsubdiv = 16

    # algoim quadrature rule degree for moment integration
    algoim_degree = 8

    # test
    w = moment_fitting_weights(ϕ, e; phase=1, npoints=npoints, algoim_degree=algoim_degree, algoim_nsubdiv=algoim_nsubdiv)
    @test isapprox(sum(w), 4-π; atol=10e-13)

    w = moment_fitting_weights(ϕ, e; phase=-1, npoints=npoints, algoim_degree=algoim_degree, algoim_nsubdiv=algoim_nsubdiv)
    @test isapprox(sum(w), π; atol=10e-13)

    w = moment_fitting_weights(ϕ, e; phase=0, npoints=npoints, algoim_degree=algoim_degree, algoim_nsubdiv=algoim_nsubdiv)
    @test isapprox(sum(w), π; atol=10e-8)
end

@testset "Moment fitting on partition" begin
    # domain
    I = Interval(0.0,2.0) ⨱ Interval(0.0,2.0)

    # number of elements
    nelem = (11,7)

    # partition
    Δ = IncreasingRange(I.data[1]..., nelem[1]+1) ⨱ IncreasingRange(I.data[2]..., nelem[2]+1)

    # levelset
    ϕ(x) = sqrt(x[1]^2 + x[2]^2) - 2.0

    # number of quadrature points
    npoints = (3,3)

    # number of element subdivisions for moment integration
    algoim_nsubdiv = 16

    # algoim quadrature rule degree for moment integration
    algoim_degree = 8

    # test
    w = moment_fitting_weights(ϕ, Δ; phase=1, npoints=npoints, algoim_degree=algoim_degree, algoim_nsubdiv=algoim_nsubdiv)
    @test isapprox(sum(w), 4-π; atol=10e-14)

    w = moment_fitting_weights(ϕ, Δ; phase=-1, npoints=npoints, algoim_degree=algoim_degree, algoim_nsubdiv=algoim_nsubdiv)
    @test isapprox(sum(w), π; atol=10e-14)

    w = moment_fitting_weights(ϕ, Δ; phase=0, npoints=npoints, algoim_degree=algoim_degree, algoim_nsubdiv=algoim_nsubdiv)
    @test isapprox(sum(w), π; atol=10e-9)
end