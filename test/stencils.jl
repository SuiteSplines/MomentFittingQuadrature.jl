@testset "Stencil dimension on single element" begin
    # element
    e = Element(Interval(0.0,2.0) ⨱ Interval(0.0,2.0))

    # stencil
    x = get_stencil(e; npoints=(7,11))

    # test dimension
    @test size(x) == (7,11)
end

@testset "Stencil dimension on Cartesian product partition" begin
    # domain
    I = Interval(0.0,2.0) ⨱ Interval(0.0,2.0)

    # number of elements
    nelem = (5,3)

    # partition
    Δ = IncreasingRange(I.data[1]..., nelem[1]+1) ⨱ IncreasingRange(I.data[2]..., nelem[2]+1)

    # number of quadrature points
    npoints = (7,11)

    # stencil
    x = get_stencil(Δ; npoints=npoints)

    # test dimension
    @test size(x) == npoints .* nelem
end

