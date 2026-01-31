using MomentFittingQuadrature, SortedSequences, StaticArrays

@testset "Single Bernstein basis function" begin
    # test single basis function evaluations
    @test bernstein(0, 0, 0.3) ≈ 1.0
    @test bernstein(1, 0, 0.0) ≈ 1.0
    @test bernstein(1, 1, 1.0) ≈ 1.0
    @test bernstein(2, 1, 0.5) ≈ 0.5

    # test edge cases
    @test bernstein(5, 0, 0.0) ≈ 1.0
    @test bernstein(5, 5, 1.0) ≈ 1.0
    @test bernstein(5, 1, 0.0) ≈ 0.0
    @test bernstein(5, 4, 1.0) ≈ 0.0

    # test partition of unity
    n = 5
    x = 0.42
    s = sum(bernstein(n, k, x) for k in 0:n)
    @test s ≈ 1.0

    # test non-negativity
    for x in IncreasingRange(0.0, 1.0, 42)
        for k in 0:n
            @test bernstein(n, k, x) ≥ 0
        end
    end
end


@testset "Mapped, single Bernstein basis function" begin
    # mapped interval [0,1] -> I
    I = Interval(-4.0, 2.0)

    # midpoint
    x = (I.a + I.b) / 2

    # test midpoint value
    n, k = 5, 3
    @test bernstein(I, n, k, x) ≈ bernstein(n, k, 0.5)

    # test endpoint values
    @test bernstein(I, n, 0, I.a) ≈ 1.0
    @test bernstein(I, n, n, I.b) ≈ 1.0

    # test partition of unity
    x = 0.42
    s = sum(bernstein(I, n, k, x) for k in 0:n)
    @test s ≈ 1.0
end


@testset "Vector of Bernstein basis functions at x" begin
    I = Interval(0.0, 1.0)
    n = 5
    x = 0.42

    # test dimension
    A = bernstein(I, n, x)
    @test size(A) == (n+1,)

    # test values
    @test all(A[k+1] ≈ bernstein(I, n, k, x) for k in 0:n)

    # test partition of unity
    @test sum(A) ≈ 1.0
end

@testset "Matrix of Bernstein basis functions at partition" begin
    I = Interval(-4.0, 2.0)
    n = 5
    x = IncreasingRange(I..., 10)

    # test dimension
    A = bernstein(I, n, x)
    @test size(A) == (n+1, length(x))

    # test each evaluation
    for l in eachindex(x)
        for k in 0:n
            @test A[k+1, l] ≈ bernstein(I, n, k, x[l])
        end
    end

    # test partition of unity columnwise
    @test all(sum(A[:, l]) ≈ 1.0 for l in 1:length(x))
end

@testset "Single tensor-product Berstein basis function at x" begin
    I = Interval(1.0, 2.0) ⨱ Interval(3.0, 4.0)
    n = (7,5)
    x = [1.7; 3.3]

    # compute univariate basis functions
    A = map(bernstein, I.data, n, x)

    # test tensor-product evaluation
    bernstein(I,n,x) ≈ kron(A[2], A[1])

    # test partition of unity
    @test sum(bernstein(I,n,x)) ≈ 1
end

@testset "Interpolation matrix of tensor-product Berstein" begin
    I = Interval(1.0, 2.0) ⨱ Interval(3.0, 4.0)
    n = (7,5)
    x = IncreasingRange(I.data[1]..., 7) ⨱ IncreasingRange(I.data[2]..., 11)

    # compute Berstein interpolation matrix
    A = bernstein(I,n,x)

    # test against single tensor-product Berstein evaluation
    @test all(A[:,k] ≈ bernstein(I,n,SVector(x[k])) for k in 1:length(x))

    # test partition of unity per column
    @test all(sum(A[:,k]) ≈ 1 for k in 1:length(x))
end