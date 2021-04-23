module GRN

using LinearAlgebra
using CSV
using DataFrames
using Distributions
using LightGraphs: SimpleDiGraph, add_edge!
using Missings
using SnowyOwl
using SparseArrays

import GLM: fit, predict, coef, stderror, loglikelihood

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
    add_moments!,

    # preprocess
    make_graph,

    # utils
    make_mapping,
    truncat_adjl!,
    truncat_gene2num!,

    # regression
    LinearRegression,
    std,
    design_matrix,
    predict,
    fit, fit!,
    residual,
    likelihood,

    # mixture
    MixtureRegression,
    hard_split,
    probabilistic_split,
    maximize_likelihood!,
    update_expectation!,
    fit!,
    fit,
    loglikelihood,
    aic,
    
    # validation
    cross_val_score,
    grid_search


include("distributions.jl")
include("layers.jl")
include("data.jl")
include("preprocess.jl")
include("utils.jl")

include("regression.jl")
include("mixture.jl")
include("validation.jl")

end # module
