### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 4e812556-727d-11eb-0b57-d11d77cfb457
begin
	using CSV, DataFrames
	using JLD2
	
	const PROJECT_PATH = dirname(@__DIR__)
	const DATA_PATH = joinpath(PROJECT_PATH, "data", "fantom5_cat")
end;

# ╔═╡ 34530a14-7282-11eb-3d8c-7d45681ffeb9
html"""
<style>
	main {
		max-width: 100%;
	}
</style>
"""

# ╔═╡ 2f515c6e-727d-11eb-2a63-c9eef8f6276f
md"# Prepare datasets"

# ╔═╡ 5a974474-727d-11eb-2f05-5f8ee74a5530
md"## Load data"

# ╔═╡ 60c600a6-727d-11eb-3a50-e3fb64b7f805
begin
	meta = CSV.File(joinpath(DATA_PATH, "time_course_meta.csv")) |> DataFrame
	first(meta, 10)
end

# ╔═╡ 94a94632-727d-11eb-00f2-d77cdee0961d
begin
	raw = CSV.File(joinpath(DATA_PATH, "time_course_raw.csv")) |> DataFrame
	first(raw, 10)
end

# ╔═╡ aa2f1c3e-727d-11eb-0121-d1bce4ca3135
md"## Select gene set"

# ╔═╡ ad3ef19a-727d-11eb-308b-3bb3efa13326
@load "../results/gene_set.jld2" gene_set

# ╔═╡ ad033402-727d-11eb-1876-49b8cae07687
nrow(raw)

# ╔═╡ acba4f62-727d-11eb-102f-89041207e066
filter!(:gene_id => x -> x in gene_set, raw)

# ╔═╡ ac7a9f34-727d-11eb-3256-7dacf5a3695c
md"## Select experiments"

# ╔═╡ ac241c86-727d-11eb-2980-5f14587130c7
exprs = ["Embryoid body to melanocyte", "ES to cardiomyocyte", "Macrophage infection",
	"MCF7 response", "MSC to adipocyte"]

# ╔═╡ abc439ea-727d-11eb-28c8-3700ed5cbad9
filter!(:experiment => x -> x in exprs, meta)

# ╔═╡ ab2de190-727d-11eb-0a62-d595a76fa19e
selected_meta = unique(meta[:, [:cell_type, :condition, :replicate]])

# ╔═╡ 415c9530-7280-11eb-286e-e9caed228bdf
nrow(selected_meta)

# ╔═╡ 436ad6c8-7280-11eb-0cb7-8fe45794cfa7
md"## Select samples"

# ╔═╡ 68625d98-7280-11eb-2f57-c140f01497d9
["H9 Embryoid body cells", "HES3-GFP Embryonic Stem cells", "Monocyte-derived macrophages",
    "MCF7 breast cancer cell line", "adipose-derived mesenchymal stem cells"]

# ╔═╡ 70785886-7280-11eb-2154-950c1e68183e
submeta = select(filter(:cell_type => (==)("adipose-derived mesenchymal stem cells"), meta), [:condition, :time, :replicate, :library_id])

# ╔═╡ df760948-7280-11eb-1d09-39b43edf3876
filter!(:replicate => (==)("biol_rep1"), submeta)

# ╔═╡ 05088f96-7281-11eb-0423-f1bf4238e441
begin
	datasets = Dict{String,Tuple}()
	exper = "adipose-derived mesenchymal stem cells"
	for rep in ["biol_rep1", "biol_rep2", "biol_rep3"]
		libs = submeta[(submeta.condition .== "adipogenic induction") .& (submeta.replicate .== rep), :library_id]
		times = submeta[submeta.replicate .== rep, :time]
		libs = append!([:gene_id], Symbol.(libs))
		datasets[exper*","*rep] = (raw[:, libs], times)
	end
	datasets
end

# ╔═╡ 23dd9240-7281-11eb-365b-e188fae35749
md"## Parsing time"

# ╔═╡ 2d0fdaa8-7281-11eb-255c-09b9bd4ee672
Set(meta.time)

# ╔═╡ 33c034ba-7281-11eb-3fe7-b15d6f635bec
keys(datasets)

# ╔═╡ 3d7acfa8-7281-11eb-39bb-01b039f72a8f
begin
	function parse_time(t::String)
		m = match(r"day(\d{2})|(\d{2})hr(?:(\d{2})min)?", t)
		if m.captures[1] != nothing
			return Meta.parse(m.captures[1]) * 1440
		else
			x = Meta.parse(m.captures[2]) * 60
			if m.captures[3] != nothing
				x += Meta.parse(m.captures[3])
			end
			return x
		end
	end
	parse_time(t) = t
end;

# ╔═╡ 4dcb890e-7281-11eb-1d11-db0aa17de169
begin
	for (k, v) in datasets
		datasets[k] = (v[1], parse_time.(v[2]))
	end
	datasets
end

# ╔═╡ 59795010-7281-11eb-08fc-971b02d1fe08
@save "../results/datasets.jld2" datasets

# ╔═╡ Cell order:
# ╟─34530a14-7282-11eb-3d8c-7d45681ffeb9
# ╟─2f515c6e-727d-11eb-2a63-c9eef8f6276f
# ╠═4e812556-727d-11eb-0b57-d11d77cfb457
# ╟─5a974474-727d-11eb-2f05-5f8ee74a5530
# ╠═60c600a6-727d-11eb-3a50-e3fb64b7f805
# ╠═94a94632-727d-11eb-00f2-d77cdee0961d
# ╟─aa2f1c3e-727d-11eb-0121-d1bce4ca3135
# ╠═ad3ef19a-727d-11eb-308b-3bb3efa13326
# ╠═ad033402-727d-11eb-1876-49b8cae07687
# ╠═acba4f62-727d-11eb-102f-89041207e066
# ╟─ac7a9f34-727d-11eb-3256-7dacf5a3695c
# ╠═ac241c86-727d-11eb-2980-5f14587130c7
# ╠═abc439ea-727d-11eb-28c8-3700ed5cbad9
# ╠═ab2de190-727d-11eb-0a62-d595a76fa19e
# ╠═415c9530-7280-11eb-286e-e9caed228bdf
# ╟─436ad6c8-7280-11eb-0cb7-8fe45794cfa7
# ╠═68625d98-7280-11eb-2f57-c140f01497d9
# ╠═70785886-7280-11eb-2154-950c1e68183e
# ╠═df760948-7280-11eb-1d09-39b43edf3876
# ╠═05088f96-7281-11eb-0423-f1bf4238e441
# ╟─23dd9240-7281-11eb-365b-e188fae35749
# ╠═2d0fdaa8-7281-11eb-255c-09b9bd4ee672
# ╠═33c034ba-7281-11eb-3fe7-b15d6f635bec
# ╠═3d7acfa8-7281-11eb-39bb-01b039f72a8f
# ╠═4dcb890e-7281-11eb-1d11-db0aa17de169
# ╠═59795010-7281-11eb-08fc-971b02d1fe08
