using CDGRNs
using DataFrames
using UMAP
using Plots
using StatsPlots
gr()

## Load data

dir = joinpath(CDGRNs.PROJECT_PATH, "results", "bonemarrow")
fig_dir = joinpath(CDGRNs.PROJECT_PATH, "pics", "bonemarrow")
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

k = 9
tree, cell_clusters = build_tree(prof, true_reg_pairs, col_palette=:default)
extract_context!(cell_clusters, tree, k)

df = get_regulation_expr(prof, tfs, true_regulations)
embedding = umap(Array(df[:, 3:end])', 2, n_neighbors=20, min_dist=0.5)
df.UMAP1 = embedding[1, :]
df.UMAP2 = embedding[2, :]
df.k9 = cell_clusters.k9

plot_configs = (markersize=2, markerstrokewidth=0, dpi=300, size=(1000, 600), thickness_scaling=2, widen=false)

p = @df df scatter(:UMAP1, :UMAP2; group=:cell, xlabel="UMAP1", ylabel="UMAP2", color_palette=:glasbey_hv_n256,
    legend=:outertopright, markersize=2, markerstrokewidth=0,
    dpi=300, size=(1000, 600), thickness_scaling=2, widen=false)
savefig(p, joinpath(fig_dir, "UMAP", "umap-cell type.svg"))
savefig(p, joinpath(fig_dir, "UMAP", "umap-cell type.png"))

p = @df df scatter(:UMAP1, :UMAP2; group=:k9, xlabel="UMAP1", ylabel="UMAP2", color_palette=:glasbey_hv_n256,
    legend=:outerright, markersize=2, markerstrokewidth=0,
    dpi=300, size=(800, 600), thickness_scaling=2, widen=false)
savefig(p, joinpath(fig_dir, "UMAP", "umap-context.svg"))
savefig(p, joinpath(fig_dir, "UMAP", "umap-context.png"))

p = @df df scatter(:UMAP1, :UMAP2; zcolor=:time, c=:coolwarm, xlabel="UMAP1", ylabel="UMAP2",
               legend=false, cb=:outerright, plot_configs...)
savefig(p, joinpath(fig_dir, "UMAP", "umap-latent time.svg"))
savefig(p, joinpath(fig_dir, "UMAP", "umap-latent time.png"))


# Train CDGRNs over contexts

cortable = train_cdgrns(tfs, prof, true_regulations,
    cell_clusters.k9, [2, 3, 6, 7, 8],
    joinpath(dir, "cdgrn_k9"))


# Network entropy

cdgrn_stats = DataFrame(cntx=String[], V=Int[], E=Int[], entropy=Float64[])
for i in [2, 7, 6, 3]
    E = nrow(cortable[i])
    V = length(unique(vcat(cortable[i].tf, cortable[i].target)))
    cortable[i].dist = cor2dist.(cortable[i].ρ)
    entr = network_entropy(to_graph(cortable[i]))
    push!(cdgrn_stats, ("$i", V, E, entr))
end
cdgrn_stats[!, :order] = [1, 2, 3, 4]

p = @df cdgrn_stats plot(:order, [:V, :E],
    label=["Number of node" "Number of edge"],
    xlabel="Context", ylabel="Network size", legend=:topright,
    thickness_scaling=2, widen=false, dpi=300, size=(800, 600)
)
xticks!([1:4;], cdgrn_stats[!,:cntx])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "bonemarrow", "CDGRN", "network_size k9.png")
savefig(p, filepath)

p = @df cdgrn_stats plot(:order, :entropy, label="Network entropy",
    xlabel="Context", ylabel="Network entropy", legend=:topright,
    thickness_scaling=2, widen=false, dpi=300, size=(800, 600)
)
xticks!([1:4;], cdgrn_stats[!,:cntx])
filepath = joinpath(CDGRNs.PROJECT_PATH, "pics", "bonemarrow", "CDGRN", "network_entropy k9.png")
savefig(p, filepath)
