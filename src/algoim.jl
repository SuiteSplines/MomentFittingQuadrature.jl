export algoim_quad_data

IgaFormation.Element(I::CartesianProduct{Dim, NTuple{Dim,T}, NTuple{Dim,Interval{T}}}) where {T,Dim} = Element(CartesianIndex(ntuple(_->1,Dim)), I)

llcorner(e::T) where {T<:Element{2}} = SVector(e[1,1])
urcorner(e::T) where {T<:Element{2}} = SVector(e[2,2])
llcorner(e::T) where {T<:Element{3}} = SVector(e[1,1,1])
urcorner(e::T) where {T<:Element{3}} = SVector(e[2,2,2])

function algoim_quad_data(f::Function, e::Element; degree::Int=1, phase::Int=-1)
    ϕ = AlgoimCallLevelSetFunction(f, x->gradient(f,x)[1])
    fill_quad_data(ϕ, llcorner(e), urcorner(e), phase, degree)
end

function algoim_quad_data(f::Function, g::Function, e::Element; degree::Int=1, phase::Int=-1)
    ϕ = AlgoimCallLevelSetFunction(f, g)
    fill_quad_data(ϕ, llcorner(e), urcorner(e), phase, degree)
end