### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ eca4ba0e-7282-11eb-3e16-db7f1b073b14
begin
	using GRN
	using DataFrames
	using CSV
	using JLD2
	using SnowyOwl
	using Gadfly
	using Statistics, InformationMeasures
	using Distributions
	using GLM
	Gadfly.set_default_plot_size(8inch, 6inch)
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
md"# Analysis of $Conc_{TF}$ and $\alpha_{targ}$"

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

# ╔═╡ f5dd857a-7286-11eb-11af-4996b7ce3454
maximum(skipmissing(prof.var.fit_likelihood))

# ╔═╡ 3ab87044-7287-11eb-09c2-7d47db077a4c
minimum(skipmissing(prof.var.fit_likelihood))

# ╔═╡ 67d7a1ea-7287-11eb-13d6-cb8f9cd7ab4b
plot(x=collect(skipmissing(prof.var.fit_likelihood)),
	 Geom.histogram,
	 Guide.title("Likelihood distribution"),
	 Guide.xlabel("Counts"),
	 Guide.ylabel("Likelihoods"),
)

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
md"## Phase diagram"

# ╔═╡ 980807ca-72ba-11eb-0347-f744ea81b3f7
begin
	# gene_name = "Pcsk2"
	# gene_name = "Ank"
	gene_name = "Rps3"
	# gene_name = "Gng12"
	i = collect(1:nrow(vars))[vars.index .== gene_name][1]
end

# ╔═╡ ba3aff2a-7295-11eb-08ee-973a494d528e
md"#### Number of nonzero spliced RNA expression: $(count(s[i, :] .!= 0)) / $(size(s, 2))"

# ╔═╡ 0be5761e-7293-11eb-3140-a37320cf1784
plot(x=s[i, :], y=u[i, :], color=prof.obs.clusters,
	 Geom.point, Guide.title("$(vars[i, :index]) phase diagram"),
	 Guide.xlabel("Spliced RNA"), Guide.ylabel("Unspliced RNA"),
)

# ╔═╡ 4bd41860-7b27-11eb-0a5d-bb7d96f2493a
plot(x=u[i, :], color=prof.obs.clusters,
	 Geom.histogram,
	 Guide.xlabel("Unspliced RNA"),
)

# ╔═╡ bfd61bba-7b26-11eb-1850-350f284b6287
md"## Phase portrait curve"

# ╔═╡ d6365c1a-7b28-11eb-3081-69d9f6d06d21
begin
	α = vars[vars.index .== gene_name, :fit_alpha][1]
	β = vars[vars.index .== gene_name, :fit_beta][1]
	γ = vars[vars.index .== gene_name, :fit_gamma][1]
	u₀ = vars[vars.index .== gene_name, :fit_u0][1]
	s₀ = vars[vars.index .== gene_name, :fit_s0][1]
end

# ╔═╡ 4b94e676-7b32-11eb-0591-c102c1ceccef
begin
	df = DataFrame(X=s[i, :], Y=u[i, :], Cell=prof.obs.clusters)
	simulate = DataFrame(τ=collect(0:0.1:100))
	simulate.u = map(τ -> SnowyOwl.unspliced(τ, u₀, α, β), simulate.τ)
	simulate.s = map(τ -> SnowyOwl.spliced(τ, s₀, u₀, α, β, γ), simulate.τ)
end;

# ╔═╡ aaf05704-7b32-11eb-0285-81e360708dc8
plot(
	layer(simulate, x=:s, y=:u, Geom.line),
	# layer(simulate, x=:u, y=:s, Geom.line),
)

# ╔═╡ e8696496-7b31-11eb-112a-cf977d01a9fa
plot(layer(simulate, x=:s, y=:u, Geom.line),
	 layer(df, x=:X, y=:Y, color=:Cell, Geom.point),
	 Guide.title("$(vars[i, :index]) phase diagram"),
	 Guide.xlabel("Spliced RNA"), Guide.ylabel("Unspliced RNA"),
)

# ╔═╡ 479ca4b2-7285-11eb-1ae4-1f7f14720b86
md"## Analysis of $s_{tf}$ and $u_{targ}$"

# ╔═╡ 7fd240e8-72c5-11eb-25b9-59198fb66347
j = 8

# ╔═╡ 5187deca-72c5-11eb-3feb-5bf0ba1babf3
plot(x=tf_s[j, :],
	 y=u[i, :], 
	 color=prof.obs.clusters,
	 Geom.point, Guide.title("Relationship of s_tf and u_targ"),
	 Guide.xlabel("Spliced RNA of TF gene, $(tf_vars[j, :index])"),
	 Guide.ylabel("Unspliced RNA of target gene, $(vars[i, :index])"),
)

# ╔═╡ 0c155df8-72c6-11eb-2299-5bf61dd3c4cd
plot(x=tf_s[j, :], Scale.x_log2(),
	 y=u[i, :], Scale.y_log2(),
	 color=prof.obs.clusters,
	 Geom.point, Guide.title("Relationship of s_tf and u_targ"),
	 Guide.xlabel("log2 spliced RNA of TF gene, $(tf_vars[j, :index])"),
	 Guide.ylabel("log2 unspliced RNA of target gene, $(vars[i, :index])"),
)

# ╔═╡ Cell order:
# ╟─6e36a6d2-86d2-11eb-210a-b5589313a599
# ╟─4786d128-7283-11eb-3191-e3fa21d398bb
# ╟─0ed0b744-7284-11eb-0632-ff6c1608fe77
# ╠═eca4ba0e-7282-11eb-3e16-db7f1b073b14
# ╟─f3f518e2-7284-11eb-07ae-393f37711c5b
# ╠═ff7bc440-7284-11eb-2e0a-7f7e6a40f1ca
# ╠═f55f2fb2-72c3-11eb-08d2-1180295c91b6
# ╟─f75f0964-7286-11eb-398b-9334c420b20b
# ╠═f5dd857a-7286-11eb-11af-4996b7ce3454
# ╠═3ab87044-7287-11eb-09c2-7d47db077a4c
# ╠═67d7a1ea-7287-11eb-13d6-cb8f9cd7ab4b
# ╠═035e2d20-7288-11eb-32f9-43bc4b4a4744
# ╠═8415622a-728e-11eb-0b5c-eb0412e77da8
# ╠═f97f2b40-7293-11eb-3137-05e11185f7b3
# ╟─a1846552-7295-11eb-37df-95e5ebb08910
# ╠═980807ca-72ba-11eb-0347-f744ea81b3f7
# ╟─ba3aff2a-7295-11eb-08ee-973a494d528e
# ╠═0be5761e-7293-11eb-3140-a37320cf1784
# ╠═4bd41860-7b27-11eb-0a5d-bb7d96f2493a
# ╟─bfd61bba-7b26-11eb-1850-350f284b6287
# ╠═d6365c1a-7b28-11eb-3081-69d9f6d06d21
# ╠═4b94e676-7b32-11eb-0591-c102c1ceccef
# ╠═aaf05704-7b32-11eb-0285-81e360708dc8
# ╠═e8696496-7b31-11eb-112a-cf977d01a9fa
# ╟─479ca4b2-7285-11eb-1ae4-1f7f14720b86
# ╠═7fd240e8-72c5-11eb-25b9-59198fb66347
# ╠═5187deca-72c5-11eb-3feb-5bf0ba1babf3
# ╠═0c155df8-72c6-11eb-2299-5bf61dd3c4cd
