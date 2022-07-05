module Analysis

using Clustering
using DataFrames
using Distances
using HypothesisTests

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
include("regulations.jl")
include("tests.jl")

end