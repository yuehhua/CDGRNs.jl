### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 7c99d08e-86b2-11eb-3cb5-d353290a2cb5
begin
	using CDGRN
	using DataFrames
	using CSV
	using JLD2
	using Gadfly
	using Statistics, Distributions
	Gadfly.set_default_plot_size(8inch, 6inch)
end;

# ╔═╡ 4e5b3c24-86bb-11eb-3231-cb8c1c9be81c
html"""
<style>
	main {
		max-width: 90%;
	}
</style>
"""

# ╔═╡ 88edd516-86e8-11eb-03d1-d1d33cccdd35
md"# Difference between before/after imputed data"

# ╔═╡ 29787ac4-86b3-11eb-3f67-99cfcf40165e
dir = joinpath(CDGRN.PROJECT_PATH, "results")

# ╔═╡ f756eb18-86b2-11eb-2c36-171f5ff2f899
begin
	obs = CSV.read(joinpath(dir, "obs.tsv"), DataFrame)
	spliced = CSV.read(joinpath(dir, "spliced.tsv"), DataFrame)
	unspliced = CSV.read(joinpath(dir, "unspliced.tsv"), DataFrame)
	Ms = CSV.read(joinpath(dir, "Ms.tsv"), DataFrame)
	Mu = CSV.read(joinpath(dir, "Mu.tsv"), DataFrame)
end;

# ╔═╡ 427d3956-86ba-11eb-22e7-6f3ad6a30734
gene_name = "Rps3"

# ╔═╡ 0c66f8c0-86ba-11eb-0627-03192fba3af3
begin
	cluster = Vector(obs[:, :clusters])
	diff_s = Ms[:, gene_name] .- spliced[:, gene_name]
	diff_u = Mu[:, gene_name] .- unspliced[:, gene_name]
end;

# ╔═╡ 2d557b94-86bd-11eb-2dd7-d9b7b15a21f5
md"## Phase portrait"

# ╔═╡ e507104c-86bb-11eb-3e5a-655ca9ec79d3
plot(x=spliced[:, gene_name], y=unspliced[:, gene_name], color=cluster,
	 Geom.point, Guide.title("$(gene_name) phase diagram"),
	 Guide.xlabel("Spliced RNA"), Guide.ylabel("Unspliced RNA"),
)

# ╔═╡ b11ae0d2-86bc-11eb-217e-e17221a60c2c
plot(x=Ms[:, gene_name], y=Mu[:, gene_name], color=cluster,
	 Geom.point, Guide.title("$(gene_name) phase diagram"),
	 Guide.xlabel("Spliced RNA"), Guide.ylabel("Unspliced RNA"),
)

# ╔═╡ 13652baa-86bd-11eb-2bd5-13a0eb0df1db
md"## Difference between spliced and moment of spliced"

# ╔═╡ 5ad5775a-86ba-11eb-0baf-d9ce29185b34
plot(x=diff_s,
	 xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster,
	 Guide.title("$(gene_name) spliced RNA histogram"),
	 Guide.xlabel("Ms - Spliced RNA"), Guide.ylabel("Count"))

# ╔═╡ ead2a03c-86ba-11eb-3a9f-3d5c83d4dc6d
plot(x=diff_u,
	 xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster,
	 Guide.title("$(gene_name) unspliced RNA histogram"),
	 Guide.xlabel("Mu - Unspliced RNA"), Guide.ylabel("Count"))

# ╔═╡ 40ded722-86be-11eb-3ebc-57e2a096c9bb
md"### Different segments"

# ╔═╡ 44f3e472-86bb-11eb-0bf7-9f55f82f24e0
plot(x=diff_s[0 .< spliced[:, gene_name] .≤ 2],
	 xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster[0 .< spliced[:, gene_name] .≤ 2],
	 Guide.title("$(gene_name) spliced RNA in (0, 2] histogram"),
	 Guide.xlabel("Ms - Spliced RNA"), Guide.ylabel("Count"))

# ╔═╡ e8cddbf0-86bd-11eb-2c66-1f82def2aab3
plot(x=diff_s[2 .< spliced[:, gene_name] .≤ 4],
	 xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster[2 .< spliced[:, gene_name] .≤ 4],
	 Guide.title("$(gene_name) spliced RNA in (2, 4] histogram"),
	 Guide.xlabel("Ms - Spliced RNA"), Guide.ylabel("Count"))

# ╔═╡ fa58f968-86bd-11eb-1c49-f53dabb3ae48
plot(x=diff_s[4 .< spliced[:, gene_name] .≤ 6],
	 xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster[4 .< spliced[:, gene_name] .≤ 6],
	 Guide.title("$(gene_name) spliced RNA in (4, 6] histogram"),
	 Guide.xlabel("Ms - Spliced RNA"), Guide.ylabel("Count"))

# ╔═╡ 09064af8-86be-11eb-28de-af10d4de1eac
plot(x=diff_s[6 .< spliced[:, gene_name] .≤ 8],
	 xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster[6 .< spliced[:, gene_name] .≤ 8],
	 Guide.title("$(gene_name) spliced RNA in (6, 8] histogram"),
	 Guide.xlabel("Ms - Spliced RNA"), Guide.ylabel("Count"))

# ╔═╡ 2e3111da-86be-11eb-3fa6-a95ebfdc387f
plot(x=diff_s[8 .< spliced[:, gene_name] .≤ 10],
	xintercept=[0], Geom.vline(color=["black"]), Geom.histogram,
	 color=cluster[8 .< spliced[:, gene_name] .≤ 10],
	 Guide.title("$(gene_name) spliced RNA in (8, 10] histogram"),
	 Guide.xlabel("Ms - Spliced RNA"), Guide.ylabel("Count"))

# ╔═╡ Cell order:
# ╟─4e5b3c24-86bb-11eb-3231-cb8c1c9be81c
# ╟─88edd516-86e8-11eb-03d1-d1d33cccdd35
# ╠═7c99d08e-86b2-11eb-3cb5-d353290a2cb5
# ╠═29787ac4-86b3-11eb-3f67-99cfcf40165e
# ╠═f756eb18-86b2-11eb-2c36-171f5ff2f899
# ╠═427d3956-86ba-11eb-22e7-6f3ad6a30734
# ╠═0c66f8c0-86ba-11eb-0627-03192fba3af3
# ╟─2d557b94-86bd-11eb-2dd7-d9b7b15a21f5
# ╠═e507104c-86bb-11eb-3e5a-655ca9ec79d3
# ╠═b11ae0d2-86bc-11eb-217e-e17221a60c2c
# ╟─13652baa-86bd-11eb-2bd5-13a0eb0df1db
# ╠═5ad5775a-86ba-11eb-0baf-d9ce29185b34
# ╠═ead2a03c-86ba-11eb-3a9f-3d5c83d4dc6d
# ╟─40ded722-86be-11eb-3ebc-57e2a096c9bb
# ╠═44f3e472-86bb-11eb-0bf7-9f55f82f24e0
# ╠═e8cddbf0-86bd-11eb-2c66-1f82def2aab3
# ╠═fa58f968-86bd-11eb-1c49-f53dabb3ae48
# ╠═09064af8-86be-11eb-28de-af10d4de1eac
# ╠═2e3111da-86be-11eb-3fa6-a95ebfdc387f
