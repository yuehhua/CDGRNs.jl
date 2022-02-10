module CDGRN

using GaussianMixtures: GaussianMixture, size
using LinearAlgebra
using Missings
using SparseArrays
using Statistics

using Clustering: kmeans, kmedoids, hclust, fuzzy_cmeans, assignments
using CSV
using DataFrames
using Distances
using Distributions
using FileIO: load
using Graphs
using MultivariateStats
using PyCall
using SnowyOwl
using StatsBase
using Colors
using Plots

import Statistics: std, cor
import GLM: fit, predict, coef, stderror, loglikelihood, dof, nobs
import StatsBase: dof, nobs, fit!

const PROJECT_PATH = dirname(@__DIR__)

export
    # data
    load_data,
    add_unspliced_data!,
    add_velocity!,
    add_moments!,
    load_tfs,
    load_CHEA,
    filter_tfs,
    filter_tfs!,

    # preprocess
    make_graph,
    get_regulation_expr,

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
    # LinearRegression,
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
    correlation,

    # mixture
    # MixtureRegression,
    # hard_split,
    # probabilistic_split,
    # maximize_likelihood!,
    # update_expectation!,
    # fit!,
    # fit,

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
    best_model,

    # plot
    make_vis_data,
    plot_regulations,

    # cdgrn
    corr_table,
    ContextDependentGRN,
    train,
    evaluate!,
    cor2dist,
    to_graph,
    network_entropy


include("data.jl")
include("preprocess.jl")
include("utils.jl")

include("clustering.jl")
include("regression.jl")
include("mixture.jl")
include("metrics.jl")
include("validation.jl")

include("plots.jl")
include("cdgrn.jl")
include("graph.jl")

end # module
