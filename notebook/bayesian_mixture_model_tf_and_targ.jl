### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ eca4ba0e-7282-11eb-3e16-db7f1b073b14
begin
	using GRN
	using SnowyOwl
	
	using DataFrames
	using CSV
	using JLD2
	using Gadfly
	
	using Statistics
	using Distributions
	using Turing
	using MCMCChains
	
	Turing.setprogress!(false);
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
md"## Select target gene"

# ╔═╡ 980807ca-72ba-11eb-0347-f744ea81b3f7
begin
	gene_name = "Rps3"
	i = collect(1:nrow(vars))[vars.index .== gene_name][1]
end

# ╔═╡ 479ca4b2-7285-11eb-1ae4-1f7f14720b86
md"## Analysis of $s_{tf}$ and $u_{targ}$"

# ╔═╡ 252009ac-8abf-11eb-2d8a-e375e0b46ec1
md"good: 2, 4, 5, 6, 7, 8; bad: 1, 3, "

# ╔═╡ 7fd240e8-72c5-11eb-25b9-59198fb66347
j = 8

# ╔═╡ e44422fa-9fe5-461e-8c0f-aef10ebd9468
begin
	df = DataFrame(X=tf_s[j, :], Y=u[i, :], Cell=prof.obs.clusters, time=prof.obs.latent_time)
	df.logX = log1p.(df.X)
	df.logY = log1p.(df.Y)
end

# ╔═╡ 5187deca-72c5-11eb-3feb-5bf0ba1babf3
p1 = plot(df, x=:X, y=:Y, color=:Cell,
	 Geom.point, Guide.title("Relationship of s_tf and u_targ"),
	 Guide.xlabel("Spliced RNA of TF gene, $(tf_vars[j, :index])"),
	 Guide.ylabel("Unspliced RNA of target gene, $(vars[i, :index])"),
)

# ╔═╡ 0c155df8-72c6-11eb-2299-5bf61dd3c4cd
p2 = plot(df, x=:logX, y=:logY, color=:Cell,
	 Geom.point, Guide.title("Relationship of s_tf and u_targ"),
	 Guide.xlabel("log2 spliced RNA of TF gene, $(tf_vars[j, :index])"),
	 Guide.ylabel("log2 unspliced RNA of target gene, $(vars[i, :index])"),
	 # Coord.cartesian(xmin=-4.8, xmax=3.5, ymin=1.5, ymax=4.9)
)

# ╔═╡ cbb5c270-4889-47f3-a612-272e47f8d67f
md"## Model"

# ╔═╡ 588bba1c-f70a-4a19-af4e-68b0c8838eae
@model DPGLM(x, y, K) = begin
    N = size(x, 1)
    
	# k
    α = 1.0
    w ~ Dirichlet(K, α)
    
    k = Vector{Int}(undef, N)
    for i in 1:N
        k[i] ~ Categorical(w)
    end
	# k ~ DirichletMultinomial(K, fill(α, N))
	
	# θ_x for each k
	x_σ = Vector(undef, K)
	x_μ = Vector(undef, K)
	for i in 1:K
		x_σ[i] ~ InverseGamma(2, 3)
		x_μ[i] ~ Normal(0, sqrt(x_σ[i]))
	end
	
	# θ_x for each observation
	x_σs = [x_σ[k[i]] for i in 1:N]
	x_μs = [x_μ[k[i]] for i in 1:N]

	# θ_y for each k
	y_σ = Vector(undef, K)
	β₀ = Vector(undef, K)
	β = Vector(undef, K)
	for i in 1:K
		y_σ[i] ~ InverseGamma(2, 1)
		β₀[i] ~ Normal(0, sqrt(3))
		β[i] ~ Normal(0, sqrt(10))
	end
	
	# θ_y for each observation
	y_σs = [y_σ[k[i]] for i in 1:N]
	y_μs = [β₀[k[i]] + x[i] * β[k[i]] for i in 1:N]
	
	# Each observation generated independently
	x ~ MvNormal(x_μs, sqrt.(x_σs))
	y ~ MvNormal(y_μs, sqrt.(y_σs))
    # return k
end

# ╔═╡ 07c4f4e6-2218-435d-867e-0d86fdafefa0
md"## Training"

# ╔═╡ 44498fa5-4081-46c7-9818-182e0de4d450
begin
	K = 4
	dpglm_model = DPGLM(df.logX, df.logY, K);
	dpglm_sampler = Gibbs(PG(100, :k), NUTS(1000, 0.85, :x_μ, :x_σ, :y_σ, :β₀, :β))
	chains = sample(dpglm_model, dpglm_sampler, MCMCThreads(), 100, 4);
end

# ╔═╡ 29fafabd-127d-4726-bf74-3457fef9f70e
size(chains)

# ╔═╡ 022aca1c-cde2-466b-968c-0099d4da467c
md"## Visualize"

# ╔═╡ 4213dfea-2f9b-46bb-a7bf-cc6a37b20bca
describe(chains)

# ╔═╡ 969f6a99-6091-494a-b503-a9a41d2876d6
begin
	function coefs(chains, K)
		w = Vector{Float64}(undef, K)
		σ = Vector{Float64}(undef, K)
		β₀ = Vector{Float64}(undef, K)
		β = Vector{Float64}(undef, K)
		for k = 1:K
			w[k] = mean(chains["w[$k]"])
			σ[k] = mean(chains["σ[$k]"])
			β₀[k] = mean(chains["β₀[$k]"])
			β[k] = mean(chains["β[$k]"])
		end
		return w, σ, β₀, β
	end
	
# 	function predict(x, chains, K)
# 		w, σ, β₀, β = coefs(chains, K)
# 		n = size(x, 1)
# 		ys = Vector{Float64}(undef, n)

# 		return 
# 	end
	
# 	function landscape(x, chains, K)
# 		w, σ, β₀, β = coefs(chains, K)

# 		σs = [σ[k[i]] for i in 1:N]
# 		μs = [β₀[k[i]] + x[i] * β[k[i]] for i in 1:N]

# 		y ~ MvNormal(μs, sqrt.(σs))
# 		return Turing.logaddexp(
# 			log(w[1]) + logpdf(Normal(0., σ[1]), x), 
# 			log(w[2]) + logpdf(Normal(0., σ[2]), x)
# 		)
# 	end
end

# ╔═╡ dadf9156-ee8d-4640-853b-c576deac50aa
begin
	ws, σs, β₀s, βs = coefs(chains, K)
	fs = [x -> β₀s[k] + x*βs[k] for k = 1:K]
end

# ╔═╡ 3ac387b5-838f-4fcd-87da-81cebee2a4b5
p = plot(
	 layer(fs, -4, 4),
	 layer(df, x=:logX, y=:logY, color=:Cell, Geom.point),
	 Guide.xlabel("log1p spliced RNA of TF gene, $(tf_vars[j, :index])"),
	 Guide.ylabel("log1p unspliced RNA of target gene, $(vars[i, :index])"),
	 Coord.cartesian(xmin=0, xmax=2.5, ymin=1.3, ymax=3.5)
)

# ╔═╡ f19d5fbf-eb39-4596-9118-e5e85b787b6f
p |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "tf-gene model", "bayesian $(tf_vars[j, :index])-$(vars[i, :index]) log plot.svg"), 8inch, 6inch)

# ╔═╡ f4a51b90-d805-486e-bfae-ec156595d86f
@save "../results/chains.jld2" chains

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
# ╠═980807ca-72ba-11eb-0347-f744ea81b3f7
# ╟─479ca4b2-7285-11eb-1ae4-1f7f14720b86
# ╟─252009ac-8abf-11eb-2d8a-e375e0b46ec1
# ╠═7fd240e8-72c5-11eb-25b9-59198fb66347
# ╠═e44422fa-9fe5-461e-8c0f-aef10ebd9468
# ╠═5187deca-72c5-11eb-3feb-5bf0ba1babf3
# ╠═0c155df8-72c6-11eb-2299-5bf61dd3c4cd
# ╟─cbb5c270-4889-47f3-a612-272e47f8d67f
# ╠═588bba1c-f70a-4a19-af4e-68b0c8838eae
# ╟─07c4f4e6-2218-435d-867e-0d86fdafefa0
# ╠═44498fa5-4081-46c7-9818-182e0de4d450
# ╠═29fafabd-127d-4726-bf74-3457fef9f70e
# ╟─022aca1c-cde2-466b-968c-0099d4da467c
# ╠═4213dfea-2f9b-46bb-a7bf-cc6a37b20bca
# ╠═969f6a99-6091-494a-b503-a9a41d2876d6
# ╠═dadf9156-ee8d-4640-853b-c576deac50aa
# ╠═3ac387b5-838f-4fcd-87da-81cebee2a4b5
# ╠═f19d5fbf-eb39-4596-9118-e5e85b787b6f
# ╠═f4a51b90-d805-486e-bfae-ec156595d86f
