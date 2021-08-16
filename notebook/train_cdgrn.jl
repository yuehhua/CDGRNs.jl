using GRN
using DataFrames
using FileIO
using JLD2
using SnowyOwl
using Distances
using Clustering
using StatsBase

## Load data

dir = joinpath(GRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)

tfs = copy(prof)

GRN.filter_genes!(prof)
vars = prof.var
u = prof.layers[:Mu]

tf_set = GRN.load_tfs(joinpath(dir, "tf_set.jld2"))
GRN.filter_tfs!(tfs, tf_set)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]


# correlation analysis

total_results = load(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results")
nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)

# map to curated TF-target database

cor_pairs = corr_table(nonsingle_pairs)
regulations = load_CHEA("/media/yuehhua/Workbench/Study/Research/PhD/data/CHEA")
reg_pairs = GRN.make_pairset(regulations)
cor_pairs.is_regulation = GRN.query_pairset(cor_pairs, reg_pairs)
true_regulations = cor_pairs[cor_pairs.is_regulation, :]


# Use true regulation pairs for clustering cells

cell_clusters = DataFrame(cell=prof.obs.clusters, time=prof.obs.latent_time)
for res in nonsingle_pairs
    if (uppercase(res[:tf_name]), uppercase(res[:gene_name])) in reg_pairs
        colname = res[:tf_name] * "_" * res[:gene_name]
        cell_clusters[!, colname] = res[:clusters]
    end
end

data = Array(cell_clusters[:, 3:end])
D = pairwise(Hamming(), data, dims=1)
GRN.clustermap(D, cell_clusters.cell, filename="clustermap_true_regulations (total)")


# Cut tree

hc_col = hclust(D, linkage=:ward, branchorder=:optimal)
cell_clusters.k3 = cutree(hc_col; k=3)
cell_clusters.k4 = cutree(hc_col; k=4)
cell_clusters.k5 = cutree(hc_col; k=5)


# Investigate cell clusters

countmap(cell_clusters[cell_clusters.k3 .== 3, :cell])


# Train on Ngn3 high EP

# context = cell_clusters.k3 .== 3
# context = cell_clusters.k5 .== 3
# context = cell_clusters.k5 .== 2
# context = cell_clusters.k5 .== 4
# context = cell_clusters.k11 .== 7
# context = cell_clusters.k11 .== 8


cdgrn = train(CDGRN, tfs, prof, true_regulations, context)
evaluate!(cdgrn)

df = GRN.partial_corr(cdgrn, tfs, prof, context)


# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)

trainX = Array(df[:, 3:end])'
p = GRN.plot_3d_pca(trainX, df.cell; context=context)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "cell clusters k3-3.html")
savefig(p, filepath)


# PCA biplot
loading = DataFrame(genes=names(df)[3:end], PC1=model.proj[:,1], PC2=model.proj[:,2], x0=0., y0=0.)
p = Gadfly.plot(loading, x=:x0, y=:y0, xend=:PC1, yend=:PC2, label=:genes,
                Geom.label, Geom.segment(arrow=true))



# Visualize regulation


# k3: 3

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k3-3 regulation $target.html")
savefig(p, filepath)


# k5: 3

target = "Ccne2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "Atad2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)


# k5: 2

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-2 regulation $target.html")
savefig(p, filepath)


# k5: 4

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)

target = "Naaladl2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Elf5", target, spliced=true)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Cpe"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Nr3c1", target, spliced=true)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Vdr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Vdr", "E2f1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Vdr", "E2f1", target, spliced=true)
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Vdr", "E2f1", target, model=pathways[Symbol(target)])
filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)

target = "Rimbp2"
target = "Ghr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", "Pax6", target, model=pathways[Symbol(target)])



# k11: 7

context = cell_clusters.k11 .== 7
