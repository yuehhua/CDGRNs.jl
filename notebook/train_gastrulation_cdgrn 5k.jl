using CDGRN
using DataFrames
using FileIO
using JLD2
using SnowyOwl
using Distances
using Clustering
using StatsBase, Statistics
using HypothesisTests
using Gadfly
import Cairo

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "gastrulation_erythroid")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)

tfs = copy(prof)

CDGRN.filter_genes!(prof)
vars = prof.var
u = prof.layers[:Mu]

tf_set = CDGRN.load_tfs(joinpath(dir, "tf_set.jld2"))
CDGRN.filter_tfs!(tfs, tf_set)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]


# correlation analysis

total_results = load(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results")
nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)

# map to curated TF-target database

cor_pairs = corr_table(nonsingle_pairs)
regulations = load_CHEA("/media/yuehhua/Workbench/Study/Research/PhD/data/CHEA")
reg_pairs = CDGRN.make_pairset(regulations)
cor_pairs.is_regulation = CDGRN.query_pairset(cor_pairs, reg_pairs)
true_regulations = cor_pairs[cor_pairs.is_regulation, :]
true_reg_pairs = filter(x -> (uppercase(x[:tf_name]), uppercase(x[:gene_name])) in reg_pairs, nonsingle_pairs)


# Use true regulation pairs for clustering cells

## hard clustering
cell_clusters = DataFrame(cell=prof.obs.clusters, time=prof.obs.latent_time)
for res in true_reg_pairs
    colname = res[:tf_name] * "_" * res[:gene_name]
    cell_clusters[!, colname] = res[:clusters]
end

data = Array(cell_clusters[:, 3:end])
D = pairwise(Hamming(), data, dims=1)
hc_col = hclust(D, linkage=:ward, branchorder=:optimal)
CDGRN.clustermap(D, cell_clusters.cell, filename="clustermap_gastrulation_erythroid")

# Cut tree

k = 5
cell_clusters.k5 = cutree(hc_col; k=k)

# Calculate context-dependent gene regulatory network

pairs = unique(zip(true_regulations.tf, true_regulations.target))
context_cor = CDGRN.cor(tfs, prof, pairs, cell_clusters.k5)
context_pairs = collect(zip(context_cor.tf, context_cor.target))
context_ρs, selected_context = CDGRN.max_cor(context_cor, k)
context_cor[!, :context] = context_ρs


# global

global_ρs = CDGRN.global_cor(tfs, prof, context_pairs, map(x -> x.size, context_ρs))
test_df = innerjoin(context_cor, global_ρs, on=[:tf, :target], makeunique=true)
test_df.global_ρ = map(x -> x.ρ, test_df.global)
test_df.global_size = map(x -> x.size, test_df.global)
test_df.context_ρ = map(x -> x.ρ, test_df.context)
test_df.context_size = map(x -> x.size, test_df.context)
plot_df = DataFrame(
    ρ=vcat(test_df.context_ρ, test_df.global_ρ),
    condition=vcat(repeat([:context], length(test_df.context_ρ)), repeat([:global], length(test_df.global_ρ))),
)
plot_df[!,:alpha] .= 0.7

bincount = 30
Gadfly.with_theme(style(grid_line_style=:solid)) do
    p = plot(plot_df, x=:ρ, color=:condition, alpha=:alpha,
            Geom.histogram(position=:identity, bincount=bincount), Guide.xlabel("TF-target pair correlation"))
    draw(SVG("pics/gastrulation_erythroid/CDGRN/context-global histogram.svg", 6inch, 4inch), p)
    draw(PNG("pics/gastrulation_erythroid/CDGRN/context-global histogram.png", 6inch, 4inch), p)
end
MannWhitneyUTest(abs.(test_df.context_ρ), abs.(test_df.global_ρ))

## KS test
r = -1.0:0.1:1.0
cntx_cdf = ecdf(test_df.context_ρ)
cntx_cdf = DataFrame(x=r, y=cntx_cdf(r), condition=:context)
global_cdf = ecdf(test_df.global_ρ)
global_cdf = DataFrame(x=r, y=global_cdf(r), condition=:global)
cdf_df = vcat(cntx_cdf, global_cdf)
Gadfly.with_theme(style(grid_line_style=:solid)) do
    p = plot(cdf_df, x=:x, y=:y, color=:condition,
            Geom.step, Guide.xlabel("TF-target pair correlation"), Guide.yticks(ticks=[0., 0.25, 0.5, 0.75, 1.]))
    draw(SVG("pics/gastrulation_erythroid/CDGRN/context-global cdf.svg", 6inch, 4inch), p)
    draw(PNG("pics/gastrulation_erythroid/CDGRN/context-global cdf.png", 6inch, 4inch), p)
end
ApproximateTwoSampleKSTest(test_df.context_ρ, test_df.global_ρ)


# spliced only

spliced_ρs = CDGRN.spliced_cor(tfs, prof, context_pairs, cell_clusters.k5, selected_context)
test_df2 = innerjoin(context_cor, spliced_ρs, on=[:tf, :target], makeunique=true)
test_df2.spliced_ρ = map(x -> x.ρ, test_df2.spliced)
test_df2.spliced_size = map(x -> x.size, test_df2.spliced)
test_df2.context_ρ = map(x -> x.ρ, test_df2.context)
test_df2.context_size = map(x -> x.size, test_df2.context)
test_df2 = test_df2[test_df2.tf .!= test_df2.target, :]
test_df2 = test_df2[.!isnan.(test_df2.spliced_ρ), :]
plot_df2 = DataFrame(
    ρ=vcat(test_df2.context_ρ, test_df2.spliced_ρ),
    condition=vcat(repeat(["unspliced+spliced"], length(test_df2.context_ρ)), repeat(["spliced"], length(test_df2.spliced_ρ))),
)
plot_df2[!,:alpha] .= 0.7

bincount = 30
Gadfly.with_theme(style(grid_line_style=:solid)) do
    p = plot(plot_df2, x=:ρ, color=:condition, alpha=:alpha,
             Geom.histogram(position=:identity, bincount=bincount), Guide.xlabel("TF-target pair correlation"))
    draw(SVG("pics/gastrulation_erythroid/CDGRN/unspliced-spliced histogram.svg", 6inch, 4inch), p)
    draw(PNG("pics/gastrulation_erythroid/CDGRN/unspliced-spliced histogram.png", 6inch, 4inch), p)
end
MannWhitneyUTest(abs.(test_df2.context_ρ), abs.(test_df2.spliced_ρ))

## KS test
r = -1.0:0.1:1.0
cntx_cdf = ecdf(test_df2.context_ρ)
cntx_cdf = DataFrame(x=r, y=cntx_cdf(r), condition="unspliced+spliced")
spliced_cdf = ecdf(test_df2.spliced_ρ)
spliced_cdf = DataFrame(x=r, y=spliced_cdf(r), condition="spliced")
cdf_df = vcat(cntx_cdf, spliced_cdf)
Gadfly.with_theme(style(grid_line_style=:solid)) do
    p = plot(cdf_df, x=:x, y=:y, color=:condition,
            Geom.step, Guide.xlabel("TF-target pair correlation"), Guide.yticks(ticks=[0., 0.25, 0.5, 0.75, 1.]))
    draw(SVG("pics/gastrulation_erythroid/CDGRN/unspliced-spliced cdf.svg", 6inch, 4inch), p)
    draw(PNG("pics/gastrulation_erythroid/CDGRN/unspliced-spliced cdf.png", 6inch, 4inch), p)
end
ApproximateTwoSampleKSTest(test_df2.context_ρ, test_df2.spliced_ρ)


# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)

trainX = Array(df[:, 3:end])'
p = CDGRN.plot_3d_pca(trainX, cell_clusters.k5)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "bonemarrow", "CDGRN", "cell clusters k5.html")
savefig(p, filepath)

# Investigate cell clusters

countmap(cell_clusters[cell_clusters.k3 .== 3, :cell])


# Train on Ngn3 high EP

context = cell_clusters.k5 .== 3
# context = cell_clusters.k5 .== 2
# context = cell_clusters.k5 .== 4


cdgrn = train(CDGRN, tfs, prof, true_regulations, context)
evaluate!(cdgrn)

context_cor = CDGRN.cor(cdgrn, tfs, prof, context)
context_cor = context_cor[.!isnan.(context_cor.ρ), :]
context_cor.regulation_type = context_cor.ρ .> 0
context_cor.regulation_strength = abs.(context_cor.ρ)
CSV.write(joinpath(dir, "cdgrn_k5", "cdgrn_k5-2.csv"), context_cor)

## Effective gene set
c = context_cor[context_cor.regulation_strength .≥ 0.3, :]
gene_set = DataFrame(gene=unique!(vcat(c.tf, c.target)))
CSV.write(joinpath(dir, "cdgrn_k5", "cdgrn_k5-2_gene_set.csv"), gene_set)


# sg = CDGRN.build_graph(c)
# CDGRN.network_entropy(sg)

# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)

trainX = Array(df[:, 3:end])'
p = CDGRN.plot_3d_pca(trainX, df.cell, context)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "cell clusters k3-3.html")
savefig(p, filepath)


# PCA biplot
loading = DataFrame(genes=names(df)[3:end], PC1=model.proj[:,1], PC2=model.proj[:,2], x0=0., y0=0.)
p = Gadfly.plot(loading, x=:x0, y=:y0, xend=:PC1, yend=:PC2, label=:genes,
                Geom.label, Geom.segment(arrow=true))



# Visualize regulation



# k5: 3

target = "Ccne2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "Atad2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)


# k5: 2

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-2 regulation $target.html")
savefig(p, filepath)


# k5: 4

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)

target = "Naaladl2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Elf5", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Cpe"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Nr3c1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Vdr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Vdr", "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)

target = "Rimbp2"
target = "Ghr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", "Pax6", target, model=cdgrn.models[Symbol(target)])
