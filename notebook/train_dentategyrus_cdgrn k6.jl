using CDGRN
using SnowyOwl
using Gadfly

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "dentategyrus")
prof = load_profile(dir)
tf_set = CDGRN.load_tfs(joinpath(dir, "tf_set.jld2"))
tfs = select_genes!(copy(prof), tf_set)

select_high_likelihood!(prof)
vars = prof.var
u = prof.layers[:Mu]

select_high_likelihood!(tfs, min_likelihood=-Inf)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]


filename = joinpath(dir, "GMM-model-selection-result.jld2")
cor_pairs, nonsingle_pairs = regulation_correlation(filename)
true_regulations, true_reg_pairs = remove_spurious_pairs(cor_pairs, nonsingle_pairs)

k = 6
tree, cell_clusters = build_tree(prof, true_reg_pairs, save="clustermap_dentategyrus")
extract_context!(cell_clusters, tree, k)

# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)
p = CDGRN.plot_3d_pca(Array(df[:, 3:end])', cell_clusters.k6)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "cell clusters k6.html")
# p = CDGRN.plot_3d_pca(Array(df[:, 3:end])', df.cell, context)
# filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "dentategyrus", "CDGRN", "cell clusters k6-2.html")
savefig(p, filepath)


# Train CDGRNs over contexts

cortable = train_cdgrns(tfs, prof, true_regulations,
    cell_clusters.k6, [1, 2, 3],
    joinpath(dir, "cdgrn_k6"))


# Network entropy

# edge count
nrow(cortable[1])
#node count
length(unique(vcat(cortable[1].tf, cortable[1].target)))

cortable[1].dist = cor2dist.(cortable[1].ρ)
cortable[2].dist = cor2dist.(cortable[2].ρ)
cortable[3].dist = cor2dist.(cortable[3].ρ)
g = to_graph(cortable[1])
network_entropy(g)



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
