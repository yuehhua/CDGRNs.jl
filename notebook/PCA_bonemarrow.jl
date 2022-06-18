using CDGRN
using DataFrames
using CSV
using JLD2
using FileIO
using MultivariateStats
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
pc = CDGRN.pca_transform(Array(df[:, 3:end])', dims=5)
df.PC1 = pc[:, 1]
df.PC2 = pc[:, 2]
df.PC3 = pc[:, 3]
df.PC4 = pc[:, 4]
df.PC5 = pc[:, 5]

plot_configs = (markersize=2, markerstrokewidth=0, dpi=300)

@df df scatter(:PC1, :PC2; group=:cell, xlabel="PC1", ylabel="PC2", color_palette=:glasbey_hv_n256, legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc12-cell type.svg"))
savefig(joinpath(fig_dir, "PCA", "pc12-cell type.png"))

@df df scatter(:PC1, :PC2; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2",
               legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc12-latent time.svg"))
savefig(joinpath(fig_dir, "PCA", "pc12-latent time.png"))

@df df scatter(:PC1, :PC3; group=:cell, xlabel="PC1", ylabel="PC3", color_palette=:glasbey_hv_n256, legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc13-cell type.svg"))
savefig(joinpath(fig_dir, "PCA", "pc13-cell type.png"))

@df df scatter(:PC1, :PC3; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC3",
               legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc13-latent time.svg"))
savefig(joinpath(fig_dir, "PCA", "pc13-latent time.png"))

@df df scatter(:PC2, :PC3; group=:cell, xlabel="PC2", ylabel="PC3", color_palette=:glasbey_hv_n256, legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc23-cell type.svg"))
savefig(joinpath(fig_dir, "PCA", "pc23-cell type.png"))

@df df scatter(:PC2, :PC3; zcolor=:time, c=:coolwarm, xlabel="PC2", ylabel="PC3",
               legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc23-latent time.svg"))
savefig(joinpath(fig_dir, "PCA", "pc23-latent time.png"))


# 3D scatter plot
plotly()
plot_configs = (markersize=1, markerstrokewidth=0)

@df df scatter(:PC1, :PC2, :PC3; group=:cell, xlabel="PC1", ylabel="PC2", zlabel="PC3", color_palette=:glasbey_hv_n256, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc123-cell type.html"))

@df df scatter(:PC1, :PC2, :PC3; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2", zlabel="PC3",
               legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc123-latent time.html"))
