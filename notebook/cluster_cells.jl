### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ eca4ba0e-7282-11eb-3e16-db7f1b073b14
begin
	using GRN
	using DataFrames
	using Distances
	using CSV
	using JLD2
	using SnowyOwl
	using Clustering
	using Statistics
	using Plots
	using StatsPlots
	gr()
	default(size = (800, 600))
end;

# ╔═╡ 6e36a6d2-86d2-11eb-210a-b5589313a599
html"""
<style>
	main {
		max-width: 90%;
	}
</style>
"""

# ╔═╡ 4786d128-7283-11eb-3191-e3fa21d398bb
md"# PCA of $Conc_{TF}$ and $\alpha_{targ}$"

# ╔═╡ 0ed0b744-7284-11eb-0632-ff6c1608fe77
md"relationship of $s_{tf}$ and $u_{targ}$"

# ╔═╡ f3f518e2-7284-11eb-07ae-393f37711c5b
md"## Load data"

# ╔═╡ ff7bc440-7284-11eb-2e0a-7f7e6a40f1ca
begin
	dir = joinpath(GRN.PROJECT_PATH, "results")
	prof = load_data(dir)
	add_unspliced_data!(prof, dir)
	add_velocity!(prof, dir)
	add_moments!(prof, dir)
end

# ╔═╡ f55f2fb2-72c3-11eb-08d2-1180295c91b6
@load "../results/tf_set.jld2" tf_set

# ╔═╡ f75f0964-7286-11eb-398b-9334c420b20b
md"## Select TF and target genes"

# ╔═╡ 035e2d20-7288-11eb-32f9-43bc4b4a4744
count(skipmissing(prof.var.fit_likelihood) .≥ 0.1)

# ╔═╡ 8415622a-728e-11eb-0b5c-eb0412e77da8
begin
	select_genes(x) = !ismissing(x) && x .≥ 0.1
	vars = filter(:fit_likelihood => select_genes, prof.var)
	data = prof.data[select_genes.(prof.var.fit_likelihood), :]
	u = prof.layers[:Mu][select_genes.(prof.var.fit_likelihood), :]
	vᵤ = prof.layers[:velocity_u][select_genes.(prof.var.fit_likelihood), :]
	s = prof.layers[:Ms][select_genes.(prof.var.fit_likelihood), :]
	vₛ = prof.layers[:velocity][select_genes.(prof.var.fit_likelihood), :]
	
	sort(vars, :fit_likelihood, rev=true)
end

# ╔═╡ f97f2b40-7293-11eb-3137-05e11185f7b3
begin
	select_tfs(x) = uppercase(x) in tf_set
	tf_vars = filter(:index => select_tfs, prof.var)
	tf_data = prof.data[select_tfs.(prof.var.index), :]
	tf_u = prof.layers[:Mu][select_tfs.(prof.var.index), :]
	tf_vᵤ = prof.layers[:velocity_u][select_tfs.(prof.var.index), :]
	tf_s = prof.layers[:Ms][select_tfs.(prof.var.index), :]
	tf_vₛ = prof.layers[:velocity][select_tfs.(prof.var.index), :]
	
	tf_data = tf_data[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_u = tf_u[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_vᵤ = tf_vᵤ[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_s = tf_s[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_vₛ = tf_vₛ[.!(ismissing.(tf_vars.fit_likelihood)), :]
	filter!(:fit_likelihood => x -> !ismissing(x), tf_vars)
end

# ╔═╡ a1846552-7295-11eb-37df-95e5ebb08910
md"## Select TF-target gene pairs"

# ╔═╡ d7efd693-0210-4521-bf0e-52fccc7c4c34
begin
	top_perc = 0.001
	@load "../results/model-selection-result.jld2" all_pairs
end

# ╔═╡ e33b42b7-528e-45e7-8f83-6add435c81eb
begin
	all_pairs.mean_score = Statistics.mean.(all_pairs.scores)
	sort!(all_pairs, :mean_score)
	pattern_pairs = all_pairs[all_pairs.best_k .!= 1, :]
	top_pairs = pattern_pairs[1:ceil(Int, top_perc * nrow(all_pairs)), :]
	tf_list = top_pairs.tf_name
	gene_list = top_pairs.gene_name
	top_pairs.tf_idx = [findfirst(tf_vars.index .== x) for x in tf_list]
	top_pairs.gene_idx = [findfirst(vars.index .== x) for x in gene_list]
end

# ╔═╡ e44422fa-9fe5-461e-8c0f-aef10ebd9468
begin
	df = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
	for k in 1:nrow(top_pairs)
		tf = top_pairs.tf_name[k]
		gene = top_pairs.gene_name[k]
		i = top_pairs.gene_idx[k]
		j = top_pairs.tf_idx[k]
		df[:, Symbol(tf * "_s")] = log1p.(tf_s[j, :])
		df[:, Symbol(gene * "_u")] = log1p.(u[i, :])
		best_k = top_pairs.best_k[k]
		logX = log1p.(tf_s[j, :])
		logY = log1p.(u[i, :])
		model = fit(MixtureRegression{best_k}, logX, logY)
		df[:, Symbol("clst_$k")] = string.(model.clusters)
	end
	df
end

# ╔═╡ d66265e1-2b43-473d-b751-796caec241f0
begin
	clst = [Symbol("clst_$i") for i = 1:12]
	unique(df[:, clst])
end

# ╔═╡ 9ddb0879-d858-47e1-a2b9-b985091b5d75
begin
	cell_cluster = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
	for col in clst
		cell_cluster[:, col] = parse.(Int, df[:, col])
	end
end

# ╔═╡ 943be3ff-4ff8-4af1-b330-e4dc1d200082
begin
	X = Array(cell_cluster[:, clst])
	D = pairwise(Hamming(), X, dims=1)
end;

# ╔═╡ 79958651-09d4-4a73-9252-084caa356ed8
hc_col = hclust(D, linkage=:ward, branchorder=:optimal)

# ╔═╡ 8212b7f3-8af4-4c24-a59f-8db5de921074
begin
	l = grid(2, 1, heights=[0.2,0.8])
	
	p = plot(
		layout = l,
		plot(hc_col, xticks=false, yticks=false),
		plot(
			D[hc_col.order,hc_col.order],
			st=:heatmap,
			colorbar=false,
			c=:YlGnBu_9,
		)
	)
end

# ╔═╡ b7ffa801-ac3c-47e2-9c8a-357ef6df0e40
savefig(p, "../pics/clustering/meta-clustering_transcription_states.svg")

# ╔═╡ a8c14754-be55-4fc0-b629-f0562db4e27e
plot(
	df.Cell[hc_col.order],
	st=:heatmap,
	colorbar=false,
)

# ╔═╡ 539878c3-d003-4b90-8abb-00e0ebb747bd
CSV.write("../results/clusters.tsv", cell_cluster)

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
# ╠═e44422fa-9fe5-461e-8c0f-aef10ebd9468
# ╠═d66265e1-2b43-473d-b751-796caec241f0
# ╠═9ddb0879-d858-47e1-a2b9-b985091b5d75
# ╠═943be3ff-4ff8-4af1-b330-e4dc1d200082
# ╠═79958651-09d4-4a73-9252-084caa356ed8
# ╠═8212b7f3-8af4-4c24-a59f-8db5de921074
# ╠═b7ffa801-ac3c-47e2-9c8a-357ef6df0e40
# ╠═a8c14754-be55-4fc0-b629-f0562db4e27e
# ╠═539878c3-d003-4b90-8abb-00e0ebb747bd
