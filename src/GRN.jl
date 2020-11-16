module GRN

using LinearAlgebra
using Flux
using Flux: glorot_uniform
using Flux: @functor
using GeometricFlux
using CSV
using DataFrames
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
    add_velocity!


include("layers.jl")
include("data.jl")


end # module
