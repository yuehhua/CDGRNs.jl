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

dir = joinpath(CDGRN.PROJECT_PATH, "results", "gastrulation_erythroid")
fig_dir = joinpath(CDGRN.PROJECT_PATH, "pics", "gastrulation_erythroid")
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

df = get_regulation_expr(prof, tfs, true_regulations, labels=:celltype)
trainX = Array(df[:, 3:end])'

model = fit(PCA, trainX; maxoutdim=5)
pc = MultivariateStats.transform(model, trainX)'
df.PC1 = pc[:, 1]
df.PC2 = pc[:, 2]
df.PC3 = pc[:, 3]
df.PC4 = pc[:, 4]
df.PC5 = pc[:, 5]

plot_configs = (markersize=2, markerstrokewidth=0, dpi=300)

@df df scatter(:PC1, :PC2; group=:cell, xlabel="PC1", ylabel="PC2", legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc12-cell type.svg"))
savefig(joinpath(fig_dir, "PCA", "pc12-cell type.png"))

@df df scatter(:PC1, :PC2; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2", legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc12-latent time.svg"))
savefig(joinpath(fig_dir, "PCA", "pc12-latent time.png"))

@df df scatter(:PC1, :PC3; group=:cell, xlabel="PC1", ylabel="PC3", legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc13-cell type.svg"))
savefig(joinpath(fig_dir, "PCA", "pc13-cell type.png"))

@df df scatter(:PC1, :PC3; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC3", legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc13-latent time.svg"))
savefig(joinpath(fig_dir, "PCA", "pc13-latent time.png"))

@df df scatter(:PC1, :PC4; group=:cell, xlabel="PC1", ylabel="PC4", legend=:outertopright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc14-cell type.svg"))
savefig(joinpath(fig_dir, "PCA", "pc14-cell type.png"))

@df df scatter(:PC1, :PC4; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC4", legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc14-latent time.svg"))
savefig(joinpath(fig_dir, "PCA", "pc14-latent time.png"))


# 3D scatter plot
plotly()
plot_configs = (markersize=1, markerstrokewidth=0)

@df df scatter(:PC1, :PC2, :PC3; group=:cell, xlabel="PC1", ylabel="PC2", zlabel="PC3", plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc123-cell type.html"))

@df df scatter(:PC1, :PC2, :PC3; zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2", zlabel="PC3", plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc123-latent time.html"))
