module GRN

using LinearAlgebra
using Flux
using Flux: glorot_uniform
using Flux: @functor
using GeometricFlux
using GraphSignals: AbstractFeaturedGraph, FeaturedGraph, has_graph, adjacency_list
using CSV
using DataFrames
using Distributions
using LightGraphs: SimpleDiGraph, add_edge!
using Missings
using SnowyOwl
using SparseArrays

const PROJECT_PATH = dirname(@__DIR__)

export
    # distributions
    DistributionTransformation,
    fit,
    transform,
    transform!,
    fit_transform,

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
    make_mapping,
    truncat_adjl!,
    truncat_gene2num!


include("distributions.jl")
include("layers.jl")
include("data.jl")
include("preprocess.jl")
include("utils.jl")

end # module
