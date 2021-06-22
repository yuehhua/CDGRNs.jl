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
	using MultivariateStats
	using Statistics
	using Plots
	using StatsPlots
	plotly()
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
	clst = [Symbol("clst_$i") for i = 1:2]
	unique(df[:, clst])
end

# ╔═╡ 88b8c3d5-4bdc-4097-83eb-739d850c39d5
function encode(df)
	unique(df)
	res = Int[]
	mapping = []
	for r = eachrow(df)
		key = Tuple(r)
		if key in mapping
			val = findfirst(x->x .== key, mapping)
		else
			push!(mapping, key)
			val = length(mapping)
		end
		push!(res, val)
	end
	return res
end

# ╔═╡ dcea6708-26d7-4175-aa9c-83ae1d1b8f8f
encode(df[:, clst])

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
# ╠═88b8c3d5-4bdc-4097-83eb-739d850c39d5
# ╠═dcea6708-26d7-4175-aa9c-83ae1d1b8f8f
