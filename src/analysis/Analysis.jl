module Analysis

using Statistics

using Clustering
using DataFrames
using Distances
using GaussianMixtureRegressions
using GLM
using Graphs
using HypothesisTests
using MultivariateStats
using SimpleWeightedGraphs
using SnowyOwl
using StatsBase

export 
    # clustering
    clustering,
    gmm_clustering,
    kmeans_clustering,
    assign_clusters,

    # regulations
    make_vis_data,

    # tree
    build_tree,
    extract_context!,

    # tests
    context_correlation,
    test_pmf,
    test_cdf

include("clustering.jl")
include("tree.jl")
include("regulations.jl")
include("visualization.jl")
include("correlations.jl")
include("tests.jl")
include("model.jl")
include("entropy.jl")
include("io.jl")

end