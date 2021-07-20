module GRN

using GaussianMixtures: GaussianMixture, size
using LinearAlgebra
using Missings
using SparseArrays
using Statistics

using Clustering: kmeans, kmedoids, fuzzy_cmeans, assignments
using CSV
using DataFrames
using Distances
using Distributions
using FileIO: load
using LightGraphs: SimpleDiGraph, add_edge!
using SnowyOwl
using StatsBase

import Statistics: std
import GLM: fit, predict, coef, stderror, loglikelihood, dof, nobs
import StatsBase: dof, nobs, fit!

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

    # clustering
    clustering,
    gmm_clustering,
    kmeans_clustering,
    assign_clusters,

    # regression
    LinearRegression,
    AbstractGMR,
    NullGMR,
    GMR,
    coef,
    std,
    stderror,
    ncoef,
    nobs,
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

    # metrics
    loglikelihood,
    membership,
    aic,
    bic,
    
    # validation
    validate_score,
    grid_search,
    best_result,
    best_score,
    best_model


include("distributions.jl")
include("layers.jl")
include("data.jl")
include("preprocess.jl")
include("utils.jl")

include("clustering.jl")
include("regression.jl")
include("mixture.jl")
include("metrics.jl")
include("validation.jl")

end # module
