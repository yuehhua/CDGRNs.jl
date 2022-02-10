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
# import Cairo

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "dentategyrus")
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

function calculate_correlation()

end

total_results = load(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results")
nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)

# map to curated TF-target database

cor_pairs = corr_table(nonsingle_pairs)
cor_pairs.is_tf = cor_pairs.tf .== cor_pairs.target
# function map_to_database!()

# end
regulations = load_CHEA("/media/yuehhua/Workbench/Study/Research/PhD/data/CHEA")
reg_pairs = CDGRN.make_pairset(regulations)
cor_pairs.is_regulation = CDGRN.query_pairset(cor_pairs, reg_pairs)
true_regulations = cor_pairs[cor_pairs.is_regulation .& .!cor_pairs.is_tf, :]
true_reg_pairs = filter(x -> (uppercase(x[:tf_name]), uppercase(x[:gene_name])) in reg_pairs, nonsingle_pairs)


# Use true regulation pairs for clustering cells

function cluster_regulation_pattern()
    
end
cell_clusters = DataFrame(cell=prof.obs.clusters, time=prof.obs.latent_time)
for res in true_reg_pairs
    colname = res[:tf_name] * "_" * res[:gene_name]
    cell_clusters[!, colname] = res[:clusters]
end

data = Array(cell_clusters[:, 3:end])
D = pairwise(Hamming(), data, dims=1)
hc_col = hclust(D, linkage=:ward, branchorder=:optimal)
CDGRN.clustermap(D, cell_clusters.cell, filename="clustermap_dentategyrus")

# Cut tree

k = 6
cell_clusters.k6 = cutree(hc_col; k=k)

# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)
trainX = Array(df[:, 3:end])'

p = CDGRN.plot_3d_pca(trainX, cell_clusters.k6)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "cell clusters k6.html")
savefig(p, filepath)

# Visualize context
p = CDGRN.plot_3d_pca(trainX, df.cell, context)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "cell clusters k6-2.html")
savefig(p, filepath)

# Investigate cell clusters

countmap(cell_clusters[cell_clusters.k6 .== 3, :cell])


# neuroblast to mature

cortable = Dict()
for c in [1, 2, 3]
    context = cell_clusters.k6 .== c
    cdgrn = train(ContextDependentGRN, tfs, prof, true_regulations, context)
    evaluate!(cdgrn)

    context_cor = CDGRN.cor(cdgrn, tfs, prof, context)
    context_cor = context_cor[.!isnan.(context_cor.ρ), :]
    context_cor.reg_type = context_cor.ρ .> 0
    context_cor.reg_stng = abs.(context_cor.ρ)
    CSV.write(joinpath(dir, "cdgrn_k6", "cdgrn_k6-$(c).csv"), context_cor)

    ## Effective gene set
    correlated = context_cor[context_cor.reg_stng .≥ 0.3, :]
    gene_set = DataFrame(gene=unique!(vcat(correlated.tf, correlated.target)))
    CSV.write(joinpath(dir, "cdgrn_k6", "cdgrn_k6-$(c)_gene_set.csv"), gene_set)

    ## add gene expression
    context_cor.tf_s = [mean(vec(get_gene_expr(tfs, String(g), :Ms))[context]) for g in context_cor.tf]
    context_cor.target_u = [mean(vec(get_gene_expr(prof, String(g), :Mu))[context]) for g in context_cor.target]
    cortable[c] = context_cor
end

## join tables
joined = outerjoin(cortable[1], cortable[2], on=[:tf, :target], renamecols="_cntx1"=>"_cntx2")
joined = outerjoin(joined, cortable[3], on=[:tf, :target], renamecols=""=>"_cntx3")
replace!(joined.ρ_cntx1, missing=>0)
replace!(joined.reg_stng_cntx1, missing=>0)
replace!(joined.tf_s_cntx1, missing=>0)
replace!(joined.target_u_cntx1, missing=>0)
replace!(joined.ρ_cntx2, missing=>0)
replace!(joined.reg_stng_cntx2, missing=>0)
replace!(joined.tf_s_cntx2, missing=>0)
replace!(joined.target_u_cntx2, missing=>0)
replace!(joined.ρ_cntx3, missing=>0)
replace!(joined.reg_stng_cntx3, missing=>0)
replace!(joined.tf_s_cntx3, missing=>0)
replace!(joined.target_u_cntx3, missing=>0)
CSV.write(joinpath(dir, "cdgrn_k6", "cdgrn_k6_all.csv"), joined)

# edge count
nrow(cortable[1])
#node count
length(unique(vcat(cortable[1].tf, cortable[1].target)))

cortable[1].dist = cor2dist.(cortable[1].ρ)
cortable[2].dist = cor2dist.(cortable[2].ρ)
cortable[3].dist = cor2dist.(cortable[3].ρ)
g = to_graph(cortable[1])
network_entropy(g)


# PCA biplot
loading = DataFrame(genes=names(df)[3:end], PC1=model.proj[:,1], PC2=model.proj[:,2], x0=0., y0=0.)
p = Gadfly.plot(loading, x=:x0, y=:y0, xend=:PC1, yend=:PC2, label=:genes,
                Geom.label, Geom.segment(arrow=true))



# Visualize regulation

# k6: 2

target = "Dlc1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Tcf3", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-2 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-3 regulation $target (global).html")
savefig(p, filepath)


# k6: 3

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-2 regulation $target.html")
savefig(p, filepath)


# k6: 1

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-4 regulation $target.html")
savefig(p, filepath)

target = "Naaladl2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-4 regulation $target.html")
savefig(p, filepath)


target = "Cpe"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "k6-4 regulation $target.html")
savefig(p, filepath)
