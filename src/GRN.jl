module GRN

using LinearAlgebra
using Flux
using Flux: glorot_uniform
using Flux: @functor
using GeometricFlux
using GraphSignals: AbstractFeaturedGraph, FeaturedGraph, has_graph, adjacency_list
using CSV
using DataFrames
using LightGraphs: SimpleDiGraph, add_edge!
using Missings
using SnowyOwl
using SparseArrays

export
    # layers
    Concentration,
    GeneRegulatory,
    update_batch_edge,
    apply_batch_message,
    propagate,

    # data
    load_data,
    add_unspliced_data!,
    add_velocity!,

    # preprocess
    make_graph,

    # utils
    make_mapping


include("layers.jl")
include("data.jl")
include("preprocess.jl")
include("utils.jl")

end # module
