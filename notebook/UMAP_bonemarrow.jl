using CDGRN
using DataFrames
using CSV
using JLD2
using FileIO
using UMAP
using Plots
using StatsPlots
gr()

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "bonemarrow")
fig_dir = joinpath(CDGRN.PROJECT_PATH, "pics", "bonemarrow")
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

df = get_regulation_expr(prof, tfs, true_regulations)
embedding = umap(Array(df[:, 3:end])', 2, n_neighbors=20, min_dist=0.5)
df.UMAP1 = embedding[1, :]
df.UMAP2 = embedding[2, :]

plot_configs = (markersize=2, markerstrokewidth=0, dpi=300)

@df df scatter(:UMAP1, :UMAP2; group=:cell, xlabel="UMAP1", ylabel="UMAP2", color_palette=:glasbey_hv_n256, legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "UMAP", "umap-cell type.svg"))
savefig(joinpath(fig_dir, "UMAP", "umap-cell type.png"))

@df df scatter(:UMAP1, :UMAP2; zcolor=:time, c=:coolwarm, xlabel="UMAP1", ylabel="UMAP2",
               legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "UMAP", "umap-latent time.svg"))
savefig(joinpath(fig_dir, "UMAP", "umap-latent time.png"))

CSV.write(joinpath(dir, "umap.tsv"), df, delim='\t')
