using CDGRN
using DataFrames
using Distances
using CSV
using JLD2
using SnowyOwl
using MultivariateStats
using Statistics
using Plots
using StatsPlots
gr()
default(size = (800, 600))

# PCA of $Conc_{TF}$ and $\alpha_{targ}$

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "pancreas")
fig_dir = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model")
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

df = get_regulation_expr(prof, tfs, true_regulations)
pc = CDGRN.pca_transform(Array(df[:, 3:end])', dims=5)
df.PC1 = pc[:, 1]
df.PC2 = pc[:, 2]
df.PC3 = pc[:, 3]
df.PC4 = pc[:, 4]
df.PC5 = pc[:, 5]

plot_configs = (markersize=2, markerstrokewidth=0, dpi=300)

@df df scatter(:PC1, :PC2, group=:cell, xlabel="PC1", ylabel="PC2", legend=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc12-cell type.svg"))

@df df scatter(:PC1, :PC2, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2", legend=false, cb=:outerright)
savefig(joinpath(fig_dir, "PCA", "pc12-latent time.svg"))

@df df scatter(:PC1, :PC3; group=:cell, xlabel="PC1", ylabel="PC3", legend_position=:bottomright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc13-cell type.svg"))

@df df scatter(:PC1, :PC3, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC3", legend=false, cb=:outerright)
savefig(joinpath(fig_dir, "PCA", "pc13-latent time.svg"))

@df df scatter(:PC1, :PC4; group=:cell, xlabel="PC1", ylabel="PC4", legend_position=:bottomright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc14-cell type.svg"))

@df df scatter(:PC1, :PC4, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC4",
			   legend=false, cb=:outerright, plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc14-latent time.svg"))


# 3D scatter plot
plotly()

@df df scatter(:PC1, :PC2, :PC3; group=:cell, xlabel="PC1", ylabel="PC2", zlabel="PC3", plot_configs...)
@df df scatter(:PC1, :PC4, :PC3; group=:cell, xlabel="PC1", ylabel="PC4", zlabel="PC3", plot_configs...)
@df df scatter(:PC1, :PC2, :PC4; group=:cell, xlabel="PC1", ylabel="PC2", zlabel="PC4", plot_configs...)
savefig(joinpath(fig_dir, "PCA", "pc123-cell type.html"))

@df df scatter(:PC1, :PC2, :PC3, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2", zlabel="PC3")
savefig(joinpath(fig_dir, "PCA", "pc123-latent time.html"))

# ╔═╡ 8e20100c-ee55-4f3b-8f28-0a6c0014c925
md"## Relationship between regulatory and cell trajectory"

# ╔═╡ 93588254-49e2-4335-a364-3d47b507abb6
begin
	gene_name = "Rps3"
	i = collect(1:nrow(vars))[vars.index .== gene_name][1]
	j = 8
	df2 = DataFrame(X=tf_s[j, :], Y=u[i, :], Cell=prof.obs.clusters, time=prof.obs.latent_time)
	df2.logX = log1p.(df2.X)
	df2.logY = log1p.(df2.Y)
	df2.Ratio = df2.logX ./ df2.logY
	df2
end

# ╔═╡ bf1a3224-422a-4a7c-903c-ab2947588765
begin
	@df df2 scatter(:logX, :logY, zcolor=:Ratio,
					c=:thermometer,
					markerstrokewidth=0,
		 			xlabel="log2 spliced RNA of TF gene, $(tf_vars[j, :index])",
		 			ylabel="log2 unspliced RNA of target gene, $(vars[i, :index])",
	)
end

# ╔═╡ 4735a5a4-50e9-46c1-b7f0-cf06f4caf4a8
savefig(joinpath(fig_dir, "PCA", "TF-gene regulation.svg"))

# ╔═╡ d37b33fa-03b6-43dc-8939-816fa67342a7
@df df scatter(:PC1, :PC2, zcolor=df2.Ratio,
			   c=:thermometer,
			   markerstrokewidth=0,
			   xlabel="PC1", ylabel="PC2"
)

# ╔═╡ 67ac0c95-24e3-4ef0-a109-6cb55380bd91
savefig(joinpath(fig_dir, "PCA", "pc12-regulation.svg"))

# ╔═╡ 1f5db7fa-b713-43cc-a0d6-364b55fd4a6b
@df df scatter(:PC1, :PC2, :PC3, zcolor=df2.Ratio,
			   c=:thermometer,
			   markerstrokewidth=0,
			   xlabel="PC1", ylabel="PC2", zlabel="PC3"
)

# ╔═╡ 9e2114cc-78ef-4b9f-b003-64ba396aaabb
savefig(joinpath(fig_dir, "PCA", "pc123-regulation.svg"))

# ╔═╡ a7005d26-98b3-438f-a106-321163e071c5
savefig(joinpath(fig_dir, "PCA", "pc123-regulation.html"))
