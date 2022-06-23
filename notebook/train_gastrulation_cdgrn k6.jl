using CDGRNs
using SnowyOwl
using DataFrames
using Plots
using StatsPlots
gr()

## Load data

dir = joinpath(CDGRNs.PROJECT_PATH, "results", "gastrulation_erythroid")
prof = load_profile(dir)
tf_set = CDGRNs.load_tfs(joinpath(dir, "tf_set.jld2"))
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
filename = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "clustering", "clustermap_gastrulation_erythroid.png")
tree, cell_clusters = build_tree(prof, true_reg_pairs, celltype=:celltype, col_palette=:default, save=filename)
extract_context!(cell_clusters, tree, k)

# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations, labels=:celltype)
p = CDGRNs.plot_3d_pca(Array(df[:, 3:end])', cell_clusters.k6)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "cell clusters k6.html")

context = cell_clusters.k6 .== 1
p = CDGRNs.plot_2d_pca(Array(df[:, 3:end])', df.cell, context)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "cell clusters k6-1.png")
savefig(p, filepath)


# Train CDGRNs over contexts

cortable = train_cdgrns(tfs, prof, true_regulations,
    cell_clusters.k6, collect(1:k),
    joinpath(dir, "cdgrn_k6"))


# Network entropy

cdgrn_stats = DataFrame(cntx=String[], V=Int[], E=Int[], entropy=Float64[])
for i in 1:k
    E = nrow(cortable[i])
    V = length(unique(vcat(cortable[i].tf, cortable[i].target)))
    cortable[i].dist = cor2dist.(cortable[i].ρ)
    entr = network_entropy(to_graph(cortable[i]))
    push!(cdgrn_stats, ("context $i", V, E, entr))
end
cdgrn_stats[!, :order] = [1, 2, 3, 4, 5, 6]

p = @df cdgrn_stats plot(:order, [:V, :E],
    label=["Number of node" "Number of edge"],
    xlabel="Context", ylabel="Network size"
)
xticks!([1:6;], cdgrn_stats[!,:cntx])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "network_size k6.png")
savefig(p, filepath)

p = @df cdgrn_stats plot(:order, :entropy,
    xlabel="Context", ylabel="Network entropy", legend=false)
xticks!([1:6;], cdgrn_stats[!,:cntx])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "network_entropy k6.png")
savefig(p, filepath)




# k3: 3

target = "Ccne2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target (global).html")
savefig(p, filepath)

target = "Atad2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target (global).html")
savefig(p, filepath)

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-3 regulation $target (global).html")
savefig(p, filepath)


# k3: 2

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-2 regulation $target.html")
savefig(p, filepath)


# k3: 4

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target.html")
savefig(p, filepath)

target = "Naaladl2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Elf5", target, spliced=true)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target (global).html")
savefig(p, filepath)


target = "Cpe"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Nr3c1", target, spliced=true)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target (global).html")
savefig(p, filepath)


target = "Vdr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Vdr", "E2f1", target, spliced=true)
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "gastrulation_erythroid", "CDGRN", "k3-4 regulation $target (global).html")
savefig(p, filepath)

target = "Rimbp2"
target = "Ghr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", "Pax6", target, model=cdgrn.models[Symbol(target)])
