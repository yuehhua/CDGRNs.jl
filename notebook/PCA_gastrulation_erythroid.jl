using CDGRN
using DataFrames
using CSV
using JLD2
using FileIO
using MultivariateStats
using Plots
using StatsPlots
plotly()

# PCA of $Conc_{TF}$ and $\alpha_{targ}$

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

model = fit(PCA, trainX; maxoutdim=3)
pc = MultivariateStats.transform(model, trainX)'
df.PC1 = pc[:, 1]
df.PC2 = pc[:, 2]
df.PC3 = pc[:, 3]

@df df scatter(:PC1, :PC2, group=:cell, xlabel="PC1", ylabel="PC2", legend=:outertopright)
savefig(joinpath(fig_dir, "PCA", "pc12-cell type.html"))

@df df scatter(:PC1, :PC2, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC2")
savefig(joinpath(fig_dir, "PCA", "pc12-latent time.html"))

@df df scatter(:PC1, :PC3, group=:cell, xlabel="PC1", ylabel="PC3", legend=:outertopright)
savefig(joinpath(fig_dir, "PCA", "pc13-cell type.html"))

@df df scatter(:PC1, :PC3, zcolor=:time, c=:coolwarm, xlabel="PC1", ylabel="PC3")
savefig(joinpath(fig_dir, "PCA", "pc13-latent time.html"))

@df df scatter(:PC1, :PC2, :PC3, group=:cell, xlabel="PC1", ylabel="PC2", zlabel="PC3")
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

# ╔═╡ Cell order:
# ╟─6e36a6d2-86d2-11eb-210a-b5589313a599
# ╟─4786d128-7283-11eb-3191-e3fa21d398bb
# ╟─0ed0b744-7284-11eb-0632-ff6c1608fe77
# ╠═eca4ba0e-7282-11eb-3e16-db7f1b073b14
# ╟─f3f518e2-7284-11eb-07ae-393f37711c5b
# ╠═ff7bc440-7284-11eb-2e0a-7f7e6a40f1ca
# ╠═f55f2fb2-72c3-11eb-08d2-1180295c91b6
# ╟─f75f0964-7286-11eb-398b-9334c420b20b
# ╠═035e2d20-7288-11eb-32f9-43bc4b4a4744
# ╠═8415622a-728e-11eb-0b5c-eb0412e77da8
# ╠═f97f2b40-7293-11eb-3137-05e11185f7b3
# ╟─a1846552-7295-11eb-37df-95e5ebb08910
# ╠═d7efd693-0210-4521-bf0e-52fccc7c4c34
# ╠═e33b42b7-528e-45e7-8f83-6add435c81eb
# ╠═980807ca-72ba-11eb-0347-f744ea81b3f7
# ╠═e44422fa-9fe5-461e-8c0f-aef10ebd9468
# ╠═aba85b62-6bcf-4566-9608-e2f930829ff5
# ╠═7388d7be-d9a1-42da-baf1-26c1a0e597f0
# ╠═ad13416f-6d2e-4672-9567-535ea584a682
# ╠═0014145d-9296-4f9e-93f6-c3fc3350d9d2
# ╠═c23b896a-b6b5-4889-9465-99e82136a6d4
# ╠═d43bd1a7-6e01-45b3-9feb-717c0fa3f9c6
# ╠═e1383e10-c1d2-4ab0-b356-8c55c251aa7d
# ╠═da9eb528-db8c-49dd-9f03-5258ec80413c
# ╠═9f37e00f-83b7-46e9-bb80-5cd86aa17624
# ╠═fa8c0932-0b48-412e-b4e0-7a91513eebb7
# ╠═3c16a1e5-d9b3-42d5-a37c-c4012689147d
# ╠═cc7a3093-5c9f-486d-8053-a3a8c103846b
# ╠═5445d95c-306d-4a80-82d8-3b66769effc4
# ╟─8e20100c-ee55-4f3b-8f28-0a6c0014c925
# ╠═93588254-49e2-4335-a364-3d47b507abb6
# ╠═bf1a3224-422a-4a7c-903c-ab2947588765
# ╠═4735a5a4-50e9-46c1-b7f0-cf06f4caf4a8
# ╠═d37b33fa-03b6-43dc-8939-816fa67342a7
# ╠═67ac0c95-24e3-4ef0-a109-6cb55380bd91
# ╠═1f5db7fa-b713-43cc-a0d6-364b55fd4a6b
# ╠═9e2114cc-78ef-4b9f-b003-64ba396aaabb
# ╠═a7005d26-98b3-438f-a106-321163e071c5
