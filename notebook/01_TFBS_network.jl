### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ cba12a86-7273-11eb-1796-2d30229604dd
begin
	using CSV, DataFrames
	using LightGraphs
	using JLD2
	const PROJECT_PATH = dirname(@__DIR__)
	const DATA_PATH = joinpath(PROJECT_PATH, "data", "fantom5_cat")
end;

# ╔═╡ b8f858f0-7273-11eb-07f8-7f3cf2a7ba71
md"# Create TFBS network"

# ╔═╡ e496f084-7273-11eb-276a-7111b6d133ed
md"## Load data"

# ╔═╡ f52577ec-7273-11eb-009e-cd7d67d77904
begin
	gene_tf = CSV.File(joinpath(DATA_PATH, "gene_tf_map.csv")) |> DataFrame
	first(gene_tf, 10)
end

# ╔═╡ 469bfc0c-7274-11eb-278f-afd6954314a5
begin
	tf = CSV.File(joinpath(DATA_PATH, "tf_list.csv")) |> DataFrame
	first(tf, 10)
end

# ╔═╡ 76adea40-7274-11eb-2c02-63498aa34e5e
begin
	total_tf_set = Set(tf.hgnc_id)
	length(total_tf_set)
end

# ╔═╡ 8a48f702-7274-11eb-0694-01c256d0472f
md"## Make HGNC mapping to gene id"

# ╔═╡ 9aa77e5c-7274-11eb-1899-05bc4bc1d51b
begin
	hgnc2id = Dict(tf[i, :hgnc_id] => tf[i, :gene_id] for i = 1:nrow(tf))
end;

# ╔═╡ ae1aa19e-7274-11eb-1995-4f81620e5691
begin
	# not all TF has hgnc id
	filter!(:tf_hgnc_id => x -> x in total_tf_set, gene_tf)
end

# ╔═╡ dd25f326-7274-11eb-0c83-b1a84ecd530f
begin
	gene_tf.tf_id = map(x -> hgnc2id[x], gene_tf.tf_hgnc_id)
	first(gene_tf, 10)
end

# ╔═╡ 4b6da60c-7279-11eb-3db9-715829c36440
begin
	hgnc_gene_set = Set(gene_tf.gene_id)
	length(hgnc_gene_set)
end

# ╔═╡ 5651dafe-7279-11eb-14ee-f117785d30b8
begin
	hgnc_tf_set = Set(gene_tf.tf_id)
	length(hgnc_tf_set)
end

# ╔═╡ 97d5e81c-7279-11eb-0ec9-95fa2ef7656a
md"## Check if gene set includes all tf set"

# ╔═╡ a24cd350-7279-11eb-2de3-b1c4a1695102
hgnc_tf_set ⊆ hgnc_gene_set

# ╔═╡ 3d0cce10-727b-11eb-3df9-13c2a93ea999
md"## Make sorted gene list"

# ╔═╡ 481e997a-727b-11eb-0fbc-65c30cbda6c4
hgnc_gene_tf_set = sort([x for x = hgnc_tf_set ∪ hgnc_gene_set])

# ╔═╡ a9e4e164-727b-11eb-2355-73215d760d96
gene2num = Dict(x => i for (i, x) in enumerate(hgnc_gene_tf_set))

# ╔═╡ b8cebed4-727b-11eb-352d-31e8f5764fa3
md"## Make simple directed graph"

# ╔═╡ c1d18e3a-727b-11eb-3429-9d631ac2caa8
dg = SimpleDiGraph(length(hgnc_gene_tf_set))

# ╔═╡ d1ffc92a-727b-11eb-318f-dbb0ff59bd56
begin
	for i = 1:nrow(gene_tf)
		g = gene2num[gene_tf[i, :gene_id]]
		tf = gene2num[gene_tf[i, :tf_id]]
		add_edge!(dg, tf, g)
	end
end

# ╔═╡ fdbf0dc8-727b-11eb-2e44-0f2310211635
dg

# ╔═╡ d1c1ff50-727b-11eb-2298-c9233f4b18e1
badjlist = dg.badjlist

# ╔═╡ f10c0af4-727b-11eb-23ef-fd297a93f1c5
begin
	@save "../results/tf_gene_network.jld2" dg
	@save "../results/gene_set.jld2" gene_set=hgnc_gene_tf_set
end

# ╔═╡ Cell order:
# ╟─b8f858f0-7273-11eb-07f8-7f3cf2a7ba71
# ╠═cba12a86-7273-11eb-1796-2d30229604dd
# ╟─e496f084-7273-11eb-276a-7111b6d133ed
# ╠═f52577ec-7273-11eb-009e-cd7d67d77904
# ╠═469bfc0c-7274-11eb-278f-afd6954314a5
# ╠═76adea40-7274-11eb-2c02-63498aa34e5e
# ╟─8a48f702-7274-11eb-0694-01c256d0472f
# ╠═9aa77e5c-7274-11eb-1899-05bc4bc1d51b
# ╠═ae1aa19e-7274-11eb-1995-4f81620e5691
# ╠═dd25f326-7274-11eb-0c83-b1a84ecd530f
# ╠═4b6da60c-7279-11eb-3db9-715829c36440
# ╠═5651dafe-7279-11eb-14ee-f117785d30b8
# ╟─97d5e81c-7279-11eb-0ec9-95fa2ef7656a
# ╠═a24cd350-7279-11eb-2de3-b1c4a1695102
# ╟─3d0cce10-727b-11eb-3df9-13c2a93ea999
# ╠═481e997a-727b-11eb-0fbc-65c30cbda6c4
# ╠═a9e4e164-727b-11eb-2355-73215d760d96
# ╟─b8cebed4-727b-11eb-352d-31e8f5764fa3
# ╠═c1d18e3a-727b-11eb-3429-9d631ac2caa8
# ╠═d1ffc92a-727b-11eb-318f-dbb0ff59bd56
# ╠═fdbf0dc8-727b-11eb-2e44-0f2310211635
# ╠═d1c1ff50-727b-11eb-2298-c9233f4b18e1
# ╠═f10c0af4-727b-11eb-23ef-fd297a93f1c5
