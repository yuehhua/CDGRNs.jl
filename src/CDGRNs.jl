module CDGRNs

using CSV
using DataFrames
using GaussianMixtureRegressions
using Graphs
using SimpleWeightedGraphs
using SnowyOwl

const PROJECT_PATH = dirname(@__DIR__)

export
    # model
    corr_table,
    ContextDependentGRN,
    train,
    evaluate!,
    cor2dist,
    to_graph,
    network_entropy,
    train_cdgrns


include("preprocess/Preprocess.jl")
include("plots/Plots.jl")
include("analysis/Analysis.jl")

include("utils.jl")
include("model.jl")
include("io.jl")

end # module
