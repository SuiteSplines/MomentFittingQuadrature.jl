export weights_patch_reshape, moment_fitting_weights

weights_patch_permute_dims(::Val{2}) = (1,3,2,4)
weights_patch_permute_dims(::Val{3}) = (1,4,2,5,3,6)

function weights_patch_reshape(weights::Matrix{T}, npoints::NTuple{Dim,Int}, nelem::NTuple{Dim,Int}) where {Dim,T}
    @assert size(weights, 1) == prod(npoints) && size(weights, 2) == prod(nelem)
    blocks = reshape(weights, npoints..., nelem...)
    permuted = permutedims(blocks, weights_patch_permute_dims(Val(Dim)))
    reshape(permuted, (nelem .* npoints)...)
end

function moment_fitting_cache(; npoints::NTuple{Dim,Int}) where {Dim}
    # reference element data
    ref_interval = CartesianProduct(ntuple(_ -> Interval(0.0, 1.0), Dim)...)
    ref_element = Element(ref_interval)
    ref_stencil = get_stencil(ref_element; npoints=npoints)

    # bernstein basis evaluated at reference stencil
    A = bernstein(ref_interval, npoints .- 1, ref_stencil)

    # moment vector cache
    b = zeros(size(A,1))

    # return precomputed data
    return A, b
end

function moment_fitting_weights!(b::V, A::M, ϕ::Function, e::Element{Dim}; npoints::NTuple{Dim,Int}, phase::Int=-1, algoim_degree::Int=maximum(npoints)+1, algoim_nsubdiv::Int=0) where {Dim, V<:AbstractVector, M<:AbstractMatrix}
    # reset moment vector cache
    b .= 0

    # get element domain
    Ω = get_element_domain(e)

    # integrate over element domain
    ω = CartesianProduct(ntuple(k -> IncreasingRange(Ω.data[k]..., 2+algoim_nsubdiv), Dim)...)
    for ee in Elements(ω)
        qx, qw = algoim_quad_data(ϕ, ee; degree=algoim_degree, phase=phase)
        for k in eachindex(qx)
            b += qw[k] .* bernstein(get_element_domain(ee), npoints .- 1, qx[k])
        end
    end

    # compute moment fitting weights
    return A \ b
end

function moment_fitting_weights(ϕ::Function, e::Element{Dim}; kwargs...) where {Dim}
    # get precomputed data and weights cache
    A, b = moment_fitting_cache(; npoints=kwargs[:npoints])

    # compute moment fitting weights
    return moment_fitting_weights!(b, A, ϕ, e; kwargs...)
end

function moment_fitting_weights(ϕ::Function, Δ::CartesianProduct{Dim}; kwargs...) where {Dim}
    # get precomputed data and weights cache
    A, b = moment_fitting_cache(; npoints=kwargs[:npoints])

    # elements iterator
    elements = Elements(Δ)

    # number of elements
    nelem = length(elements)

    # weights per element cache
    w = zeros(prod(kwargs[:npoints]), nelem)    

    # loop over all elements in partition Δ and compute moment fitting weights
    for (eind,e) in enumerate(elements)
        w[:,eind] = moment_fitting_weights!(b, A, ϕ, e; kwargs...)
    end

    # return weights
    return w
end