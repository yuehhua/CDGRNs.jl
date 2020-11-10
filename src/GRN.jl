module GRN

using LinearAlgebra
using Flux
using Flux: glorot_uniform
using Flux: @functor
using GeometricFlux

export
    # layers
    Concentration,
    GeneRegulatory,
    update_batch_edge,
    apply_batch_message,
    propagate


include("layers.jl")


end # module
