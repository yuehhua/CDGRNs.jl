using GRN
using DataFrames
using Distances
using CSV
using JLD2
using SnowyOwl
using UMAP
using Statistics
using Plots
using StatsPlots
plotly()
default(size = (800, 600))

dir = joinpath(GRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)

@load "results/tf_set.jld2" tf_set

select_genes(x) = !ismissing(x) && x .≥ 0.1
vars = filter(:fit_likelihood => select_genes, prof.var)
data = prof.data[select_genes.(prof.var.fit_likelihood), :]
u = prof.layers[:Mu][select_genes.(prof.var.fit_likelihood), :]
vᵤ = prof.layers[:velocity_u][select_genes.(prof.var.fit_likelihood), :]
s = prof.layers[:Ms][select_genes.(prof.var.fit_likelihood), :]
vₛ = prof.layers[:velocity][select_genes.(prof.var.fit_likelihood), :]

sort(vars, :fit_likelihood, rev=true)

select_tfs(x) = uppercase(x) in tf_set
tf_vars = filter(:index => select_tfs, prof.var)
tf_data = prof.data[select_tfs.(prof.var.index), :]
tf_u = prof.layers[:Mu][select_tfs.(prof.var.index), :]
tf_vᵤ = prof.layers[:velocity_u][select_tfs.(prof.var.index), :]
tf_s = prof.layers[:Ms][select_tfs.(prof.var.index), :]
tf_vₛ = prof.layers[:velocity][select_tfs.(prof.var.index), :]

tf_data = tf_data[.!(ismissing.(tf_vars.fit_likelihood)), :]
tf_u = tf_u[.!(ismissing.(tf_vars.fit_likelihood)), :]
tf_vᵤ = tf_vᵤ[.!(ismissing.(tf_vars.fit_likelihood)), :]
tf_s = tf_s[.!(ismissing.(tf_vars.fit_likelihood)), :]
tf_vₛ = tf_vₛ[.!(ismissing.(tf_vars.fit_likelihood)), :]
filter!(:fit_likelihood => x -> !ismissing(x), tf_vars)

manual = false
top_perc = 0.001
manual || @load "results/model-selection-result.jld2" all_pairs
folder_name = manual ? "pca-trajectory-manual" : "pca-trajectory"

if manual
    tf_list = ["E2f1", "Nr3c1", "Pax6", "Pdx1", "Vdr"]
    gene_list = ["Auts2", "Bicc1", "Cacna1d", "Cadps", "Cdc14b", "Cdh1", "Cntrob", "Cpe", "Dcbld1", "Dnmt3a", "Ezh2", "Fam219a", "Pcsk2", "Rps3"]
else
    all_pairs.mean_score = Statistics.mean.(all_pairs.scores)
    sort!(all_pairs, :mean_score)
    pattern_pairs = all_pairs[all_pairs.best_k .!= 1, :]
    top_pairs = pattern_pairs[1:ceil(Int, top_perc * nrow(all_pairs)), :]
    tf_list = top_pairs.tf_name
    gene_list = top_pairs.gene_name
end

tf_idx = [findfirst(tf_vars.index .== x) for x in tf_list]
gene_idx = [findfirst(vars.index .== x) for x in gene_list]

df = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
for (j, tf) = zip(tf_idx, tf_list)
    for (i, gene) = zip(gene_idx, gene_list)
        df[:, Symbol(tf * "_s")] = log1p.(tf_s[j, :])
        df[:, Symbol(gene * "_u")] = log1p.(u[i, :])
    end
end

# n_components = 2
n_neighbors = 30

trainX = Array(df[:, 3:end])'
# embedding = umap(trainX, n_components, n_neighbors=n_neighbors, min_dist=0.5)
# df.UMAP1 = embedding[1, :]
# df.UMAP2 = embedding[2, :]

# @df df scatter(:UMAP1, :UMAP2, group=:Cell, xlabel="UMAP1", ylabel="UMAP2")
# savefig(joinpath(GRN.PROJECT_PATH, "pics", folder_name, "umap12-cell type.svg"))


n_components = 3
embedding = umap(trainX, n_components, n_neighbors=n_neighbors, min_dist=0.5)
df.UMAP1 = embedding[1, :]
df.UMAP2 = embedding[2, :]
df.UMAP3 = embedding[3, :]

@df df scatter(:UMAP1, :UMAP2, :UMAP3, group=:Cell, xlabel="UMAP1", ylabel="UMAP2", zlabel="UMAP3")
savefig(joinpath(GRN.PROJECT_PATH, "pics", folder_name, "umap123-cell type.html"))
