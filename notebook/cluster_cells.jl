using CDGRN
using DataFrames
using FileIO
using JLD2
# using Gadfly
using SnowyOwl
using Distances
# Gadfly.set_default_plot_size(8inch, 6inch)
# import Cairo

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)


# correlation analysis

total_results = load(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results")
nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)
correlations = map(x -> correlation(x[:model]), nonsingle_pairs)


# clustering cells against TF-gene pair components

cell_clusters = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
for res in nonsingle_pairs
    colname = res[:tf_name] * "_" * res[:gene_name]
    cell_clusters[!, colname] = res[:clusters]
end

data = Array(cell_clusters[:, 3:end])
D = pairwise(Hamming(), data, dims=1)
CDGRN.clustermap(D, cell_clusters.Cell, filename="clustermap (total)")


# map to curated TF-target database

## unit: component
cor_pairs = map(x -> (tf=x[:tf_name], target=x[:gene_name], best_k=x[:best_k], corr=correlation(x[:model])), nonsingle_pairs)
cor_pairs = flatten(DataFrame(cor_pairs), :corr)
cor_pairs.adjusted_corr = CDGRN.fisher_transform(cor_pairs.corr)

reg_pairs = Set(map(i -> (regulations[i,:tf], regulations[i,:target]), 1:nrow(regulations)))

query_pairs = map(i -> (uppercase(cor_pairs[i,:tf]), uppercase(cor_pairs[i,:target])), 1:nrow(cor_pairs))
cor_pairs.is_regulation = map(x -> x in reg_pairs, query_pairs)
true_regulations = cor_pairs[cor_pairs.is_regulation, :]


# Use true regulation pairs for clustering cells

cell_clusters = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
for res in nonsingle_pairs
    if (uppercase(res[:tf_name]), uppercase(res[:gene_name])) in reg_pairs
        colname = res[:tf_name] * "_" * res[:gene_name]
        cell_clusters[!, colname] = res[:clusters]
    end
end

data = Array(cell_clusters[:, 3:end])
D = pairwise(Hamming(), data, dims=1)
CDGRN.clustermap(D, cell_clusters.Cell, filename="clustermap_true_regulations (total)")


# Use true regulation pairs for PCA

tfs = copy(prof)

CDGRN.filter_genes!(prof)
vars = prof.var
u = prof.layers[:Mu]

tf_set = CDGRN.load_tfs(joinpath(dir, "tf_set.jld2"))
CDGRN.filter_tfs!(tfs, tf_set)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]

tf_list = unique(true_regulations.tf)
gene_list = unique(true_regulations.target)

df = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
for gene in gene_list
    df[:, Symbol(gene * "_u")] = vec(get_gene_expr(prof, gene, :Mu))
end
for tf in tf_list
    df[:, Symbol(tf * "_s")] = vec(get_gene_expr(tfs, tf, :Ms))
end

using MultivariateStats

trainX = Array(df[:, 3:end])'
model = fit(PCA, trainX; maxoutdim=3)
pc = MultivariateStats.transform(model, trainX)'
df.PC1 = pc[:, 1]
df.PC2 = pc[:, 2]
df.PC3 = pc[:, 3]

using Plots
using StatsPlots
plotly()

p = @df df scatter(:PC1, :PC2, :PC3, group=:Cell, xlabel="PC1", ylabel="PC2", zlabel="PC3")
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "PCA", "pc123-cell type_true reg.html")
savefig(p, filepath)

p = @df df scatter(:PC1, :PC2, :PC3, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2", zlabel="PC3")
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "PCA", "pc123-latent time_true reg.html")
savefig(p, filepath)


col_colors = CDGRN.generate_column_colors(cell_clusters.Cell)
col_colors = map(x -> RGB(x...), col_colors)
p = @df df scatter(:PC1, :PC2, :PC3, color=col_colors, xlabel="PC1", ylabel="PC2", zlabel="PC3")
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "clustering", "PCA with clustermap color labels.html")
savefig(p, filepath)
