### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 509f6f2e-a729-11eb-33f1-ddf731d351ec
begin
	using GRN
	using DataFrames
	using CSV
	using JLD2
	using SnowyOwl
	using GLM
	using Plots
	using StatsPlots
	plotly()
	default(size = (800, 600))
end;

# ╔═╡ b482472f-8e28-42f5-8054-51a98409935f
html"""
<style>
	main {
		max-width: 90%;
	}
</style>
"""

# ╔═╡ f9f576e8-064e-4766-83da-bb07233903d0
md"## Load data"

# ╔═╡ 48b83213-b400-4e68-9a45-62823724ed3c
begin
	dir = joinpath(GRN.PROJECT_PATH, "results")
	prof = load_data(dir)
	add_unspliced_data!(prof, dir)
	add_velocity!(prof, dir)
	add_moments!(prof, dir)
end

# ╔═╡ 03b537f3-6e84-4921-9bc5-f81634f500d3
@load "../results/tf_set.jld2" tf_set

# ╔═╡ ea80b584-d91d-4645-839c-f8430fb00eb8
md"## Select TF and target genes"

# ╔═╡ 9b1ce956-1872-4d62-a23d-b6bc9a0e6084
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

# ╔═╡ f04d0029-ce31-409b-9f1b-a666ebc05f9f
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

# ╔═╡ 68905d15-9cc2-4969-ac4b-f0e164698d8b
md"## Select TF-target gene pairs"

# ╔═╡ 2d8693ad-a015-4025-85a8-fe1041cd2c5d
begin
	tf_name = "Pdx1"
	gene_name = "Gck"
	j = findfirst(tf_vars.index .== tf_name)
	i = findfirst(vars.index .== gene_name)
end

# ╔═╡ c5d4099c-c3d8-4984-b219-251797b464fe
begin
	df = DataFrame(X=tf_s[j, :], Y=u[i, :], Cell=prof.obs.clusters, time=prof.obs.latent_time)
	df.logX = log1p.(df.X)
	df.logY = log1p.(df.Y)
end

# ╔═╡ 2a9b2a44-b2c3-4c05-8267-683846a2aa0f
@df df scatter(:X, :Y, group=:Cell,
	 xlabel="Spliced RNA of TF gene, $(tf_name)",
	 ylabel="Unspliced RNA of target gene, $(gene_name)",
)

# ╔═╡ 4000c964-18d8-4c75-90b7-09cf3acdf0d7
savefig("../pics/tf-regulation/$(tf_name)-$(gene_name) regulation.svg")

# ╔═╡ b2edd330-cd3a-4dcd-b6b4-a8317faee76d
filtered_df = filter(:Cell => x -> !(x in ["Delta", "Ngn3 high EP", "Epsilon", "Ductal", "Ngn3 low EP", "Beta"]), df)

# ╔═╡ 370fea8f-c49d-4b9e-9958-2c6c3795bc3c
@df filtered_df scatter(:X, :Y, group=:Cell,
	 xlabel="Spliced RNA of TF gene, $(tf_name)",
	 ylabel="Unspliced RNA of target gene, $(gene_name)",
)

# ╔═╡ c017b034-31bf-4ec9-96ce-6c10af204208
savefig("../pics/tf-regulation/$(tf_name)-$(gene_name) regulation-2.svg")

# ╔═╡ bb877d3b-1f26-4593-b963-bd3e115e9e35
md"## Modeling"

# ╔═╡ 6c4bcd36-9a62-471f-962d-25644f00a290
begin
	ymax = maximum(filtered_df.Y)
	p = (filtered_df.Y .+ ymax)./(2ymax)
	filtered_df.Odds = p ./ (1 .- p)
	filtered_df.logOdds = log.(filtered_df.Odds)
	filtered_df.X_inv = 1 ./ filtered_df.X
	filtered_df.Y_inv = 1 ./ filtered_df.Y
end

# ╔═╡ ddeb082f-fa24-447f-b9e9-fb9623ef3c2a
@df filtered_df scatter(:X_inv, :Y_inv, group=:Cell)

# ╔═╡ f8b8d1e0-e19c-43c0-b011-d2b7e350baf2
model = lm(@formula(Y_inv ~ X_inv), filtered_df)

# ╔═╡ cbde93d5-254e-43fd-8d47-3489a57de9d2
begin
	x = collect(0.1:0.01:7.5)
	y = predict(model, DataFrame(X_inv=x))
end

# ╔═╡ 7605af83-1901-4e2d-8f49-92f46f78171f
begin
	@df filtered_df scatter(:X_inv, :Y_inv, group=:Cell)
	plot!(x, y)
end

# ╔═╡ Cell order:
# ╟─b482472f-8e28-42f5-8054-51a98409935f
# ╠═509f6f2e-a729-11eb-33f1-ddf731d351ec
# ╟─f9f576e8-064e-4766-83da-bb07233903d0
# ╠═48b83213-b400-4e68-9a45-62823724ed3c
# ╠═03b537f3-6e84-4921-9bc5-f81634f500d3
# ╟─ea80b584-d91d-4645-839c-f8430fb00eb8
# ╠═9b1ce956-1872-4d62-a23d-b6bc9a0e6084
# ╠═f04d0029-ce31-409b-9f1b-a666ebc05f9f
# ╠═68905d15-9cc2-4969-ac4b-f0e164698d8b
# ╠═2d8693ad-a015-4025-85a8-fe1041cd2c5d
# ╠═c5d4099c-c3d8-4984-b219-251797b464fe
# ╠═2a9b2a44-b2c3-4c05-8267-683846a2aa0f
# ╠═4000c964-18d8-4c75-90b7-09cf3acdf0d7
# ╠═b2edd330-cd3a-4dcd-b6b4-a8317faee76d
# ╠═370fea8f-c49d-4b9e-9958-2c6c3795bc3c
# ╠═c017b034-31bf-4ec9-96ce-6c10af204208
# ╟─bb877d3b-1f26-4593-b963-bd3e115e9e35
# ╠═6c4bcd36-9a62-471f-962d-25644f00a290
# ╠═ddeb082f-fa24-447f-b9e9-fb9623ef3c2a
# ╠═f8b8d1e0-e19c-43c0-b011-d2b7e350baf2
# ╠═cbde93d5-254e-43fd-8d47-3489a57de9d2
# ╠═7605af83-1901-4e2d-8f49-92f46f78171f
