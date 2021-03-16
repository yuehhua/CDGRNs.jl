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
end;

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
	u = prof.layers[:unspliced][select_genes.(prof.var.fit_likelihood), :]
	vᵤ = prof.layers[:velocity_u][select_genes.(prof.var.fit_likelihood), :]
	s = prof.layers[:spliced][select_genes.(prof.var.fit_likelihood), :]
	vₛ = prof.layers[:velocity][select_genes.(prof.var.fit_likelihood), :]
	
	sort(vars, :fit_likelihood, rev=true)
end

# ╔═╡ f97f2b40-7293-11eb-3137-05e11185f7b3
begin
	select_tfs(x) = uppercase(x) in tf_set
	tf_vars = filter(:index => select_tfs, prof.var)
	tf_data = prof.data[select_tfs.(prof.var.index), :]
	tf_u = prof.layers[:unspliced][select_tfs.(prof.var.index), :]
	tf_vᵤ = prof.layers[:velocity_u][select_tfs.(prof.var.index), :]
	tf_s = prof.layers[:spliced][select_tfs.(prof.var.index), :]
	tf_vₛ = prof.layers[:velocity][select_tfs.(prof.var.index), :]
	
	tf_data = tf_data[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_u = tf_u[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_vᵤ = tf_vᵤ[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_s = tf_s[.!(ismissing.(tf_vars.fit_likelihood)), :]
	tf_vₛ = tf_vₛ[.!(ismissing.(tf_vars.fit_likelihood)), :]
	filter!(:fit_likelihood => x -> !ismissing(x), tf_vars)
end

# ╔═╡ 3acfb3b6-728b-11eb-3828-91342e1f0175
# select_prof = prof[:, skipmissing(prof.var.fit_likelihood) .≥ 0.1]

# ╔═╡ a1846552-7295-11eb-37df-95e5ebb08910
md"## Phase diagram"

# ╔═╡ 980807ca-72ba-11eb-0347-f744ea81b3f7
begin
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
	 # Coord.cartesian(xmin=0, xmax=50, ymin=0, ymax=50),
)

# ╔═╡ 589bc586-72ba-11eb-0aad-b382564a1b31
plot(x=s[i, :], y=u[i, :], color=prof.obs.clusters,
	 Geom.point, Guide.title("$(vars[i, :index]) phase diagram"),
	 Guide.xlabel("log2 spliced RNA"), Guide.ylabel("log2 unspliced RNA"),
	 Scale.x_log2(), Scale.y_log2(),
	 # Coord.cartesian(xmin=-1, xmax=6, ymin=-1, ymax=6),
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

# ╔═╡ 4167e540-72cf-11eb-18f5-4f3e17034ddb
md"## Mutual information"

# ╔═╡ 57a8fa8a-72cf-11eb-2d00-19cb2f5fd0af
get_mutual_information(u[i, :], tf_s[1, :])

# ╔═╡ e5575cd0-72cf-11eb-2184-e95277440312
get_mutual_information(u[i, :], tf_s[4, :])

# ╔═╡ ea7fdfd6-72cf-11eb-2983-053d7d4d4d08
get_mutual_information(u[i, :], tf_s[5, :])

# ╔═╡ efe9d85a-72cf-11eb-0bcf-7f06e2b17042
get_mutual_information(u[i, :], tf_s[6, :])

# ╔═╡ ffe172e2-72cf-11eb-2d41-0dca085cddec
get_mutual_information(u[i, :], tf_s[9, :])

# ╔═╡ f6c8325c-72cf-11eb-2b18-39f335eb996d
md"---"

# ╔═╡ cae5fb4c-72cf-11eb-1100-351f5303e5e1
get_mutual_information(u[i, :], tf_s[2, :])

# ╔═╡ cd0f198a-72cf-11eb-3cf4-ed34bdfe7f25
get_mutual_information(u[i, :], tf_s[3, :])

# ╔═╡ ce45a760-72cf-11eb-1f7b-d73aedcdc8ca
get_mutual_information(u[i, :], tf_s[7, :])

# ╔═╡ d2513fea-72cf-11eb-2fe5-f714f79de89c
get_mutual_information(u[i, :], tf_s[8, :])

# ╔═╡ 0f17c7f0-72d0-11eb-3616-ff5d79cbfb70
mis = [get_mutual_information(u[i, :], tf_s[k, :]) for k = 1:9]

# ╔═╡ 231babf4-72d0-11eb-3035-418885ef8fa0
median(mis)

# ╔═╡ Cell order:
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
# ╠═3acfb3b6-728b-11eb-3828-91342e1f0175
# ╟─a1846552-7295-11eb-37df-95e5ebb08910
# ╠═980807ca-72ba-11eb-0347-f744ea81b3f7
# ╟─ba3aff2a-7295-11eb-08ee-973a494d528e
# ╠═0be5761e-7293-11eb-3140-a37320cf1784
# ╠═589bc586-72ba-11eb-0aad-b382564a1b31
# ╟─479ca4b2-7285-11eb-1ae4-1f7f14720b86
# ╠═7fd240e8-72c5-11eb-25b9-59198fb66347
# ╠═5187deca-72c5-11eb-3feb-5bf0ba1babf3
# ╠═0c155df8-72c6-11eb-2299-5bf61dd3c4cd
# ╟─4167e540-72cf-11eb-18f5-4f3e17034ddb
# ╠═57a8fa8a-72cf-11eb-2d00-19cb2f5fd0af
# ╠═e5575cd0-72cf-11eb-2184-e95277440312
# ╠═ea7fdfd6-72cf-11eb-2983-053d7d4d4d08
# ╠═efe9d85a-72cf-11eb-0bcf-7f06e2b17042
# ╠═ffe172e2-72cf-11eb-2d41-0dca085cddec
# ╟─f6c8325c-72cf-11eb-2b18-39f335eb996d
# ╠═cae5fb4c-72cf-11eb-1100-351f5303e5e1
# ╠═cd0f198a-72cf-11eb-3cf4-ed34bdfe7f25
# ╠═ce45a760-72cf-11eb-1f7b-d73aedcdc8ca
# ╠═d2513fea-72cf-11eb-2fe5-f714f79de89c
# ╠═0f17c7f0-72d0-11eb-3616-ff5d79cbfb70
# ╠═231babf4-72d0-11eb-3035-418885ef8fa0
