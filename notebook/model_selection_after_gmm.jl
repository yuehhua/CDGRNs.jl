# Model selection

using GRN
using DataFrames
using CSV
using JLD2
using SnowyOwl
using Gadfly
using MLDataUtils
using Clustering
using Distances
using Statistics
Gadfly.set_default_plot_size(8inch, 6inch)

## Load data

dir = joinpath(GRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)

@load "results/tf_set.jld2" tf_set

## Select TF and target genes

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

# prepare data
gene_name = "Rps3"
i = collect(1:nrow(vars))[vars.index .== gene_name][1]
j = 8
df = DataFrame(X=tf_s[j, :], Y=u[i, :])
df.logX = log1p.(df.X)
df.logY = log1p.(df.Y)

# train
k_range = 2:5
models = [GaussianMixtures.GMM(k, train, kind=:full) for k in k_range]

# predict
lls = [GaussianMixtures.llpg(m, train) for m in models]
clusters = [vec(map(x -> x[2], argmax(ll, dims=2))) for ll in lls]

# validate
D = pairwise(Euclidean(), train, dims=1)
sil_scores = [silhouettes(clst, D) for clst in clusters]
mean_sil_scores = [mean(s) for s in sil_scores]
