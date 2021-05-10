# Model selection

using GRN
using DataFrames
using CSV
using JLD2
using SnowyOwl
using Gadfly
using MLDataUtils
using Statistics: mean
Gadfly.set_default_plot_size(8inch, 6inch)

## Load data

dir = joinpath(GRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)

@load "results/tf_set.jld2" tf_set

## Select TF and target genes

select_genes(x) = !ismissing(x) && x .≥ 0.1
vars = filter(:fit_likelihood => select_genes, prof.var)
data = prof.data[select_genes.(prof.var.fit_likelihood), :]
u = prof.layers[:Mu][select_genes.(prof.var.fit_likelihood), :]
vᵤ = prof.layers[:velocity_u][select_genes.(prof.var.fit_likelihood), :]
s = prof.layers[:Ms][select_genes.(prof.var.fit_likelihood), :]
vₛ = prof.layers[:velocity][select_genes.(prof.var.fit_likelihood), :]

sort(vars, :fit_likelihood, rev=true)

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

## All pairs combination

# Threads.@threads for i = 1:nrow(vars)
# 	for j = 1:nrow(tf_vars)
# 		df = DataFrame(X=tf_s[j, :], Y=u[i, :])
# 		df.logX = log1p.(df.X)
# 		df.logY = log1p.(df.Y)
# 		p = plot(df, x=:logX, y=:logY,
# 				Geom.point, Guide.title("Relationship of s_tf and u_targ"),
# 				Guide.xlabel("log1p spliced RNA of TF gene, $(tf_vars[j, :index])"),
# 				Guide.ylabel("log1p unspliced RNA of target gene, $(vars[i, :index])"),
# 		)
# 		p |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "all pair tf-gene", "$(tf_vars[j, :index])-$(vars[i, :index]) log plot.svg"), 10inch, 6inch)
# 	end
# end

gene_name = "Rps3"
i = collect(1:nrow(vars))[vars.index .== gene_name][1]
j = 8
df = DataFrame(X=tf_s[j, :], Y=u[i, :])
df.logX = log1p.(df.X)
df.logY = log1p.(df.Y)

k_range = 2:10
cv = 5
best_k = grid_search(MixtureRegression, df.logX, df.logY, k_range, cv=cv, verbosity=2)
