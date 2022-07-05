module CDGRNs

using DataFrames
using Graphs

const PROJECT_PATH = dirname(@__DIR__)

include("preprocess/Preprocess.jl")
include("plots/Plots.jl")
include("analysis/Analysis.jl")
include("utils.jl")

end # module
