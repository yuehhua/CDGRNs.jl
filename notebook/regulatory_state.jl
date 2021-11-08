### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ ebbc8a0a-c8b7-11eb-0a11-b568ab547e85
begin
	using CDGRN
	using DataFrames
	using Distances
	using CSV
	using JLD2
	using SnowyOwl
	using Statistics
	using Plots
	using StatsPlots
	plotly()
	default(size = (800, 600))
end

# ╔═╡ 94d45f5e-80b9-47ea-ad3b-cbf69bdbe359
html"""
<style>
	main {
		max-width: 90%;
	}
</style>
"""

# ╔═╡ 921a0994-52e6-4c82-bfe0-8fd1724a9621
md"## Load data"

# ╔═╡ 85e9c386-c64f-425f-98d8-1828d93f28ce
begin
	dir = joinpath(CDGRN.PROJECT_PATH, "results")
	prof = load_data(dir)
	add_unspliced_data!(prof, dir)
	add_velocity!(prof, dir)
	add_moments!(prof, dir)
end

# ╔═╡ 089a43a2-04f3-4528-89b6-a04c165c6829
@load "../results/tf_set.jld2" tf_set

# ╔═╡ bfffac0d-9d68-44f1-9898-e0a061494f32
@load "../results/model-selection-result.jld2" all_pairs

# ╔═╡ f8dccb50-a80e-4160-8417-91f2772b9a89
md"## Select TF and target genes"

# ╔═╡ 68b63d4a-9494-47b7-8981-2761d81ecafa
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

# ╔═╡ 64f8503b-cbfe-4de9-9d16-fa22141f2bd8
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

# ╔═╡ 3c611856-5554-4bed-b333-27ccd44c2fa8
md"## Select top 0.1% pairs"

# ╔═╡ 58c1b02f-9126-4674-9106-19885a6188e4
begin
	top_perc = 0.001
	all_pairs.mean_score = mean.(all_pairs.scores)
	sort!(all_pairs, :mean_score)
	pattern_pairs = all_pairs[all_pairs.best_k .!= 1, :]
	top_pairs = pattern_pairs[1:ceil(Int, top_perc * nrow(all_pairs)), :]
	tf_list = top_pairs.tf_name
	gene_list = top_pairs.gene_name
	tf_idx = [findfirst(tf_vars.index .== x) for x in tf_list]
	gene_idx = [findfirst(vars.index .== x) for x in gene_list]
end

# ╔═╡ 964f7cd9-f57b-427a-ae01-a4eff2744c5a
md"## Regard components in mixture model as regulatory states"

# ╔═╡ f09d2b91-1fd4-425d-a4ef-2a12b2ad0fc1
Set(top_pairs.tf_name)

# ╔═╡ c6299417-cccb-48f0-b1f9-24d6cd3288c8
Set(top_pairs.gene_name)

# ╔═╡ 60fb4d51-4477-4923-9616-9ed892b9b794
top_pairs

# ╔═╡ da96591b-d837-454c-90c8-4d3736076d5c
begin
	df = DataFrame(Cell=prof.obs.clusters, time=prof.obs.latent_time)
	for (j, tf) = zip(tf_idx, tf_list)
		for (i, gene) = zip(gene_idx, gene_list)
			df[:, Symbol(tf * "_s")] = log1p.(tf_s[j, :])
			df[:, Symbol(gene * "_u")] = log1p.(u[i, :])
		end
	end
	df
end

# ╔═╡ ca3fccde-106b-4361-b58c-887f63b8100e
begin
	pair_index = 12
	i = gene_idx[pair_index]
	j = tf_idx[pair_index]
	df2 = DataFrame(X=tf_s[j, :], Y=u[i, :], Cell=prof.obs.clusters, 	
		time=prof.obs.latent_time)
	df2.logX = log1p.(df2.X)
	df2.logY = log1p.(df2.Y)
	df2
end

# ╔═╡ 9b896a7c-a205-4d28-ad27-b8b0cdc773ac
begin
	@df df2 scatter(:logX, :logY, group=:Cell)
	savefig(joinpath(CDGRN.PROJECT_PATH, "pics", "model-selection", "top", "top-$(pair_index) $(tf_list[pair_index]) - $(gene_list[pair_index]).svg"))
end

# ╔═╡ Cell order:
# ╟─94d45f5e-80b9-47ea-ad3b-cbf69bdbe359
# ╠═ebbc8a0a-c8b7-11eb-0a11-b568ab547e85
# ╟─921a0994-52e6-4c82-bfe0-8fd1724a9621
# ╠═85e9c386-c64f-425f-98d8-1828d93f28ce
# ╠═089a43a2-04f3-4528-89b6-a04c165c6829
# ╠═bfffac0d-9d68-44f1-9898-e0a061494f32
# ╟─f8dccb50-a80e-4160-8417-91f2772b9a89
# ╠═68b63d4a-9494-47b7-8981-2761d81ecafa
# ╠═64f8503b-cbfe-4de9-9d16-fa22141f2bd8
# ╟─3c611856-5554-4bed-b333-27ccd44c2fa8
# ╠═58c1b02f-9126-4674-9106-19885a6188e4
# ╟─964f7cd9-f57b-427a-ae01-a4eff2744c5a
# ╠═f09d2b91-1fd4-425d-a4ef-2a12b2ad0fc1
# ╠═c6299417-cccb-48f0-b1f9-24d6cd3288c8
# ╠═60fb4d51-4477-4923-9616-9ed892b9b794
# ╠═da96591b-d837-454c-90c8-4d3736076d5c
# ╠═ca3fccde-106b-4361-b58c-887f63b8100e
# ╠═9b896a7c-a205-4d28-ad27-b8b0cdc773ac
