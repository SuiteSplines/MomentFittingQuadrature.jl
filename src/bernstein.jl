export bernstein

function bernstein(n::Int, k::Int, x::Real)
    return binomial(n,k) * x.^k .* (1.0 .- x).^(n-k)
end

function bernstein(I::Interval, n::Int, k::Int, y::Real)
    x = (y-I.a) / (I.b-I.a)
    return bernstein(n,k,x)
end

function bernstein(I::Interval, n::Int, x::Real)
    A = zeros(n+1)
    for i in 1:n+1
        A[i] = bernstein(I,n,i-1,x)
    end
    return A
end

function bernstein(I::Interval, n::Int, x::AbstractVector)
    A = zeros(n+1, length(x))
    for j in eachindex(x)
        for i in 1:n+1
            A[i,j] = bernstein(I,n,i-1,x[j])
        end
    end
    return A
end

function bernstein(I::CartesianProduct{Dim}, n::NTuple{Dim,Int}, x::AbstractVector) where {Dim}
    KroneckerProduct(ntuple(i -> bernstein(I.data[i], n[i], x[i]), Dim)...; reverse=true)
end

function bernstein(I::CartesianProduct{Dim}, n::NTuple{Dim,Int}, x::CartesianProduct{Dim}) where {Dim}
    KroneckerProduct(ntuple(k -> bernstein(I.data[k], n[k], x.data[k]), Dim)...; reverse=true)
end