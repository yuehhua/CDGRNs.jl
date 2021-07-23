using GRN
using DataFrames
using FileIO
using JLD2
using Gadfly
using SnowyOwl
using Distances
Gadfly.set_default_plot_size(8inch, 6inch)

## Load data

dir = joinpath(GRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)


# correlation analysis

total_results = load(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results")
nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)
correlations = map(x -> correlation(x[:model]), nonsingle_pairs)

corr_dist = GRN.fisher_transform(vcat(correlations...))
p = plot(x=corr_dist, Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist.svg")
p |> SVG(filepath, 5inch, 3inch)

abs_corr_dist = abs.(corr_dist)
count(abs_corr_dist .>= 0.5)
count(abs_corr_dist .>= 0.6)
count(abs_corr_dist .>= 0.7)
count(abs_corr_dist .>= 0.8)
count(abs_corr_dist .>= 0.9)


# clustering cells against TF-gene pair components

cell_clusters = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
for res in nonsingle_pairs
    colname = res[:tf_name] * "_" * res[:gene_name]
    cell_clusters[!, colname] = res[:clusters]
end

data = Array(cell_clusters[:, 3:end])
D = pairwise(Hamming(), data, dims=1)
GRN.clustermap(D, cell_clusters.Cell)


# Cut tree
using Clustering
using Plots
using StatsPlots
gr()

hc_col = hclust(D, linkage=:ward, branchorder=:optimal)
cell_clusters[!, :k3] = cutree(hc_col; k=3)
cell_clusters[!, :k4] = cutree(hc_col; k=4)
cell_clusters[!, :k8] = cutree(hc_col; k=8)

encode(xs) = Dict(zip(unique(xs), 1:length(unique(xs))))
mapping = encode(cell_clusters.Cell)
cells = map(x -> mapping[x], cell_clusters.Cell)
randindex(cells, cell_clusters.k3)

conting_k3 = counts(cells, cell_clusters.k3)
conting_k4 = counts(cells, cell_clusters.k4)
conting_k8 = counts(cells, cell_clusters.k8)
