### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 7f21b56a-7415-11eb-1da1-c51dc8846ce3
begin
	using CDGRN
	using DataFrames
	using CSV
	using JLD2
	using SnowyOwl
	using Gadfly
	using Statistics, Distributions
end;

# ╔═╡ 68d66740-7415-11eb-1e9e-0dc03b9e14eb
md"# Analysis of distributions of gene expression"

# ╔═╡ b3783af0-7415-11eb-3abd-87dbeb2a6eea
md"## Load data"

# ╔═╡ a7ed5576-7415-11eb-31c1-550fdf692958
begin
	dir = joinpath(CDGRN.PROJECT_PATH, "results")
	prof = load_data(dir)
	add_unspliced_data!(prof, dir)
	add_velocity!(prof, dir)
end

# ╔═╡ b831d0d8-7415-11eb-2bc4-7f8f84884315
prof.var.index

# ╔═╡ 7de93a22-74c1-11eb-136e-d59157e77a26
md"## Visualization"

# ╔═╡ f2e627e2-7415-11eb-1238-03cd2c44f694
begin
	gene_name = "Actr1b"
	xs = prof.data[prof.var.index .== gene_name, :]
	nonzero_xs = xs[xs .!= 0]
	plot(x=nonzero_xs,
		 Geom.histogram,
		 Guide.xlabel("Expression"),
	)
end

# ╔═╡ 89af577e-74c1-11eb-0269-bb10ce4be866
md"## Model"

# ╔═╡ 78708596-74c1-11eb-2f00-b5a7b891caa3
model = Distributions.fit(Gamma, collect(nonzero_xs))

# ╔═╡ 0a2f26f8-74c3-11eb-29c5-4982d5b461f0
plot(layer(x->pdf(model, x), 0, 2.5, color=[colorant"black"]),
	 layer(x=nonzero_xs, Geom.histogram(density=true)),
	 Guide.xlabel("Expression"),
)

# ╔═╡ 670bbfea-7580-11eb-27af-956a85df7a5f
md"### Correction by normal distribution"

# ╔═╡ 63a2a626-757f-11eb-1c67-43eb91ea3767
begin
	dt = CDGRN.fit(DistributionTransformation, Gamma, nonzero_xs)
	zs = CDGRN.transform(dt, nonzero_xs)
	normal = dt.dst
	plot(layer(x->pdf(normal, x), 0, 2.5, color=[colorant"black"]),
		 layer(x=zs, Geom.histogram(density=true)),
		 Guide.xlabel("Expression"),
	)
end

# ╔═╡ Cell order:
# ╟─68d66740-7415-11eb-1e9e-0dc03b9e14eb
# ╟─b3783af0-7415-11eb-3abd-87dbeb2a6eea
# ╠═7f21b56a-7415-11eb-1da1-c51dc8846ce3
# ╠═a7ed5576-7415-11eb-31c1-550fdf692958
# ╠═b831d0d8-7415-11eb-2bc4-7f8f84884315
# ╟─7de93a22-74c1-11eb-136e-d59157e77a26
# ╠═f2e627e2-7415-11eb-1238-03cd2c44f694
# ╟─89af577e-74c1-11eb-0269-bb10ce4be866
# ╠═78708596-74c1-11eb-2f00-b5a7b891caa3
# ╠═0a2f26f8-74c3-11eb-29c5-4982d5b461f0
# ╟─670bbfea-7580-11eb-27af-956a85df7a5f
# ╠═63a2a626-757f-11eb-1c67-43eb91ea3767
