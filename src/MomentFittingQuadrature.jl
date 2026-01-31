module MomentFittingQuadrature

using SortedSequences, CartesianProducts, IgaFormation, UnivariateSplines, KroneckerProducts
using StaticArrays, Algoim, Zygote

include("bernstein.jl")
include("algoim.jl")
include("stencils.jl")
include("weights.jl")

end
