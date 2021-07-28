using GRN
using DataFrames
using FileIO
using JLD2
using Gadfly
using SnowyOwl
using Distances
Gadfly.set_default_plot_size(8inch, 6inch)
import Cairo

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

# corr_dist_z3 = corr_dist[(corr_dist .>= 3) .| (corr_dist .<= -3)]
# p = plot(x=corr_dist_z3, Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
# filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist_z3.png")
# p |> PNG(filepath, 5inch, 3inch)


# map to curated TF-target database
## unit: component

cor_pairs = map(x -> (tf=x[:tf_name], target=x[:gene_name], corr=correlation(x[:model])), nonsingle_pairs)
cor_pairs = flatten(DataFrame(cor_pairs), :corr)
cor_pairs.adjusted_corr = GRN.fisher_transform(cor_pairs.corr)

reg_pairs = Set(map(i -> (regulations[i,:tf], regulations[i,:target]), 1:nrow(regulations)))

query_pairs = map(i -> (uppercase(cor_pairs[i,:tf]), uppercase(cor_pairs[i,:target])), 1:nrow(cor_pairs))
cor_pairs.is_regulation = map(x -> x in reg_pairs, query_pairs)  # 83

unique(cor_pairs[cor_pairs.is_regulation, :], [:tf, :target])  # 28

p = plot(cor_pairs, x=:adjusted_corr, color=:is_regulation,
         Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist with label.svg")
p |> SVG(filepath, 5inch, 3inch)

p = plot(x=cor_pairs[cor_pairs.is_regulation, :adjusted_corr],
         Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist (true regulation).svg")
p |> SVG(filepath, 5inch, 3inch)

## unit: pair

cor_pairs = map(x -> (tf=x[:tf_name], target=x[:gene_name], corr=correlation(x[:model])), nonsingle_pairs)
cor_pairs = DataFrame(cor_pairs)
# The maximum absolute (component) correlation as pair correlation
cor_pairs.corr = map(x -> x[argmax(abs.(x))], cor_pairs.corr)
cor_pairs.adjusted_corr = GRN.fisher_transform(cor_pairs.corr)

query_pairs = map(i -> (uppercase(cor_pairs[i,:tf]), uppercase(cor_pairs[i,:target])), 1:nrow(cor_pairs))
cor_pairs.is_regulation = map(x -> x in reg_pairs, query_pairs)

using HypothesisTests

tests = []
for z in 2:150
    xs = cor_pairs.adjusted_corr .> z
    ys = cor_pairs.is_regulation
    a = count(xs .& ys)
    b = count(.!(xs) .& ys)
    c = count(xs .& .!(ys))
    d = count(.!(xs) .& .!(ys))
    test = HypothesisTests.FisherExactTest(a, b, c, d)
    push!(tests, test)
end
pvalues = pvalue.(tests)


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
