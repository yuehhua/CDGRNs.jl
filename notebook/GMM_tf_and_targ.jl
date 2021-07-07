using Logging

using GRN
using DataFrames
using CSV
using JLD2
using SnowyOwl
using Gadfly
using Statistics
using GaussianMixtures
using Distributions
using Clustering
using Distances
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

## Select target gene

# gene_name = "Rps3"
# i = collect(1:nrow(vars))[vars.index .== gene_name][1]

## Analysis of $s_{tf}$ and $u_{targ}$

# j = 8
# df = DataFrame(X=tf_s[j, :], Y=u[i, :], Cell=prof.obs.clusters, time=prof.obs.latent_time)
# df.logX = log1p.(df.X)
# df.logY = log1p.(df.Y)
# train = hcat(df.logX, df.logY)

# l1 = layer(df, x=:X, y=:Y, color=:Cell, Geom.point)
# p1 = plot(l1,
# 	 Guide.title("Relationship of s_tf and u_targ"),
# 	 Guide.xlabel("Spliced RNA of TF gene, $(tf_vars[j, :index])"),
# 	 Guide.ylabel("Unspliced RNA of target gene, $(vars[i, :index])"),
# )

# l2 = layer(df, x=:logX, y=:logY, color=:Cell, Geom.point)
# p2 = plot(l2,
# 	 Guide.title("Relationship of s_tf and u_targ"),
# 	 Guide.xlabel("log2 spliced RNA of TF gene, $(tf_vars[j, :index])"),
# 	 Guide.ylabel("log2 unspliced RNA of target gene, $(vars[i, :index])"),
# )

# k = 4
# gmm = GaussianMixtures.GMM(k, train, kind=:full)
# posterior, ll = GaussianMixtures.gmmposterior(gmm, train)
# df.cluster = string.(assign_clusters(posterior))

# results = grid_search(df.logX, df.logY, 2:5, best_model=true, verbosity=2)
# model = results[:model]
# mix_logpdf(x,y) = logpdf(model, [x,y])

# D = pairwise(Euclidean(), X, dims=1)
# ms = mean(silhouettes(clusters, D))

k_range = 1:5
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
            data = hcat(df.logX, df.logY)

            results = grid_search(GMR, data, k_range, best_model=true, return_score=true)
            best_k = results[:best_k]
            model = results[:model]
			scores = results[:score]
			clusters = results[:clusters]
            mix_logpdf(x,y) = logpdf(model, [x,y])

            lock(splock) do
                @info "(i=$i, j=$j) $tf_name - $gene_name: best k = $best_k"
                push!(total_results, (tf_name=tf_name, gene_name=gene_name, best_k=best_k, scores=scores))
            end

            if best_k != 1
                df.clusters = string.(clusters)
                xmax = ceil(maximum(df.logX))
                xmin = floor(minimum(df.logX))
                ymax = ceil(maximum(df.logY))
                ymin = floor(minimum(df.logY))

				# log likelihood plot
				l3 = layer(df, x=:logX, y=:logY, Geom.point)
				l4 = layer(z=mix_logpdf, xmin=[xmin], xmax=[xmax], ymin=[ymin], ymax=[ymax], Geom.contour)
				coord = Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
				p3 = plot(l3, l4, coord,
					Guide.xlabel("log1p spliced RNA of TF gene, $(tf_name)"),
					Guide.ylabel("log1p unspliced RNA of target gene, $(gene_name)"),
				)
				p3 |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log likelihood plot.svg"), 10inch, 6inch)


				# cluster plot
				l5 = layer(df, x=:logX, y=:logY, color=:cluster, Geom.point)
				p4 = plot(l5, coord,
					Guide.xlabel("log1p spliced RNA of TF gene, $(tf_name)"),
					Guide.ylabel("log1p unspliced RNA of target gene, $(gene_name)"),
				)
				p4 |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log cluster plot.svg"), 10inch, 6inch)
            end
        end
    end
end

close(io)

all_pairs = DataFrame(total_results)

@save "results/GMM-model-selection-result.jld2" all_pairs
