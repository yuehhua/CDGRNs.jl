# Model selection
using Logging
using Statistics: mean

using CSV
using DataFrames
using Gadfly
using JLD2
using MLDataUtils

using GRN
using SnowyOwl

Gadfly.set_default_plot_size(8inch, 6inch)

# logger
io = open("model-selection.log", "w+")
logger = SimpleLogger(io)

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

k_range = 1:5
cv = 5
splock = Threads.SpinLock()
total_results = []
with_logger(logger) do
    Threads.@threads for i = 1:nrow(vars)
        for j = 1:nrow(tf_vars)
            tf_name = tf_vars[j, :index]
            gene_name = vars[i, :index]

            df = DataFrame(X=tf_s[j, :], Y=u[i, :])
            df.logX = log1p.(df.X)
            df.logY = log1p.(df.Y)

            results = grid_search(MixtureRegression, df.logX, df.logY, k_range,
                                  λ=0.02, cv=cv, 
                                  best_model=true, return_score=true, verbosity=2)
            best_k = results[:best_k]
            model = results[:model]
            scores = results[:score]
            ll = loglikelihood(model, average=true)

            lock(splock) do
                @info "(i=$i, j=$j) $tf_name - $gene_name: best k = $best_k with log likelihood: $ll"
                push!(total_results, (tf_name=tf_name, gene_name=gene_name, best_k=best_k, ll=ll, scores=scores))
            end

            if best_k != 1
                df.clusters = string.(model.clusters)
                fs = [x -> coef(model.models[i])'*[1, x] for i = 1:best_k]
                xmax = ceil(maximum(df.logX))
                xmin = floor(minimum(df.logX))

                p = plot(
                        layer(df, x=:logX, y=:logY, color=:clusters, Geom.point),
                        layer(fs, xmin, xmax),
                        Guide.xlabel("log1p spliced RNA of TF gene, $(tf_name)"),
                        Guide.ylabel("log1p unspliced RNA of target gene, $(gene_name)"),
                )
                p |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "model-selection", "$(tf_name)-$(gene_name) log plot.svg"), 10inch, 6inch)
            end
        end
    end
end

close(io)

all_pairs = DataFrame(total_results)

@save "results/model-selection-result.jld2" all_pairs

# gene_name = "Rps3"
# i = collect(1:nrow(vars))[vars.index .== gene_name][1]
# j = 8
# df = DataFrame(X=tf_s[j, :], Y=u[i, :])
# df.logX = log1p.(df.X)
# df.logY = log1p.(df.Y)

# best_k = grid_search(MixtureRegression, df.logX, df.logY, k_range, λ=0.02, verbosity=2)
