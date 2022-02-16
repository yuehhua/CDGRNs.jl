using CDGRN
using SnowyOwl
using Gadfly

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "pancreas")
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

k = 5
filename = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "clustering", "clustermap_pancreatic.png")
tree, cell_clusters = build_tree(prof, true_reg_pairs, col_palette=:default, save=filename)
extract_context!(cell_clusters, tree, k)


# Compare context and global

context_cor, selected_context, context_pairs = context_correlation(tfs, prof, true_regulations, cell_clusters.k5, k)
global_ρs = CDGRN.global_cor(tfs, prof, context_pairs, map(x -> x.size, context_cor.context))
test_df = innerjoin(context_cor, global_ρs, on=[:tf, :target], makeunique=true)
context_ρ = map(x -> x.ρ, test_df.context)
global_ρ = map(x -> x.ρ, test_df.global)
test_result = test_pmf(context_ρ, global_ρ, "context", "global"; plot_dir="pics/tf-gene gmm model/CDGRN", title="context-global")
test_result = test_cdf(context_ρ, global_ρ, "context", "global"; plot_dir="pics/tf-gene gmm model/CDGRN", title="context-global")


# Compare spliced only

spliced_ρs = CDGRN.spliced_cor(tfs, prof, context_pairs, cell_clusters.k5, selected_context)
test_df2 = innerjoin(context_cor, spliced_ρs, on=[:tf, :target], makeunique=true)
test_df2.spliced_ρ = map(x -> x.ρ, test_df2.spliced)
test_df2.context_ρ = map(x -> x.ρ, test_df2.context)
test_df2 = test_df2[test_df2.tf .!= test_df2.target, :]
test_df2 = test_df2[.!isnan.(test_df2.spliced_ρ), :]
test_result = test_pmf(test_df2.context_ρ, test_df2.spliced_ρ, "unspliced+spliced", "spliced"; plot_dir="pics/tf-gene gmm model/CDGRN", title="unspliced-spliced")
test_result = test_cdf(test_df2.context_ρ, test_df2.spliced_ρ, "unspliced+spliced", "spliced"; plot_dir="pics/tf-gene gmm model/CDGRN", title="unspliced-spliced")


# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)
p = CDGRN.plot_3d_pca(Array(df[:, 3:end])', cell_clusters.k5)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "cell clusters k5.html")

context = cell_clusters.k5 .== 1
p = CDGRN.plot_2d_pca(Array(df[:, 3:end])', df.cell, context)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "cell clusters k5-1.png")
savefig(p, filepath)


# Train CDGRNs over contexts

cortable = train_cdgrns(tfs, prof, true_regulations,
    cell_clusters.k5, [1, 2, 3, 4, 5],
    joinpath(dir, "cdgrn_k5"))


# Network entropy

for i in 1:5
    E = nrow(cortable[i])
    V = length(unique(vcat(cortable[i].tf, cortable[i].target)))
    cortable[i].dist = cor2dist.(cortable[i].ρ)
    entr = network_entropy(to_graph(cortable[i]))
    println("Context $i:\nnode:\t$V\nedge:\t$E\nentropy:\t$entr")
end



# k5: 3
context = cell_clusters.k5 .== 3

target = "Ccne2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "Atad2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)


# k5: 2
context = cell_clusters.k5 .== 2

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-2 regulation $target.html")
savefig(p, filepath)


# k5: 4
context = cell_clusters.k5 .== 4

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)

target = "Naaladl2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Elf5", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Cpe"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Nr3c1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Vdr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)

target = "Rimbp2"
target = "Ghr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", "Pax6", target, model=cdgrn.models[Symbol(target)])
