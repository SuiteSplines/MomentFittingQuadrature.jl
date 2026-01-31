export get_stencil

function get_stencil(I::Interval{T}; npoints::Int=1) where {T}
    quadrule = GaussRule(Legendre, npoints, I.a, I.b)
    quadrule.x
end

function get_stencil(e::Element{Dim}; npoints::NTuple{Dim,Int}) where {Dim}
    domain = get_element_domain(e)
    CartesianProduct(ntuple(k -> get_stencil(domain.data[k]; npoints=npoints[k]), Dim)...)
end

function get_stencil(Δ::IncreasingRange{T}; npoints::Int) where {T}
    # number of elements
    m = length(Δ) - 1

    # total number of quadrature points
    n = m * npoints

    # stencil cache
    x = zeros(n)

    # loop over partition
    for (k,inds) in enumerate(Iterators.partition(1:n, npoints))
        x[inds] = get_stencil(Interval(Δ[k], Δ[k+1]); npoints=npoints)
    end
    
    # return stencil
    return x
end

function get_stencil(Δ::CartesianProduct{Dim}; npoints::NTuple{Dim,Int}) where {Dim}
    #CartesianProduct(ntuple(k -> PatchRule(Δ.data[k]; npoints=npoints[k]).x[2:end-1], Dim)...)
    CartesianProduct(ntuple(k -> get_stencil(Δ.data[k]; npoints=npoints[k]), Dim)...)
end