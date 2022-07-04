module CDGRNs

using Missings

using Clustering
using CSV
using DataFrames
using Distances
# using FileIO: load
using GaussianMixtureRegressions
using Graphs
using SnowyOwl
using SimpleWeightedGraphs

const PROJECT_PATH = dirname(@__DIR__)

export
    # interface
    load_profile,
    regulation_correlation,
    remove_spurious_pairs,
    build_tree,
    extract_context!,
    context_correlation,
    test_pmf,
    test_cdf,
    train_cdgrns,

    # data
    load_data,
    add_unspliced_data!,
    add_velocity!,
    add_moments!,
    load_tfs,
    load_CHEA,
    select_high_likelihood!,
    select_genes!,

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

include("plots.jl")
include("model.jl")
include("graph.jl")

include("interface.jl")

include("plots/Plots.jl")

end # module
