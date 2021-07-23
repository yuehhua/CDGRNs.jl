using Logging

using GRN
using DataFrames
using CSV
using FileIO
using JLD2
using SnowyOwl
using Gadfly
using Statistics
using GaussianMixtures
using Distributions
using Clustering
using Distances
Gadfly.set_default_plot_size(8inch, 6inch)

## Load data

dir = joinpath(GRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)

## Select TF and target genes

select_genes(x) = !ismissing(x) && x .≥ 0.1
vars = filter(:fit_likelihood => select_genes, prof.var)
data = prof.data[select_genes.(prof.var.fit_likelihood), :]
u = prof.layers[:Mu][select_genes.(prof.var.fit_likelihood), :]

tf_set = GRN.load_tfs(joinpath(dir, "tf_set.jld2"))

select_tfs(x) = uppercase(x) in tf_set
tf_vars = filter(:index => select_tfs, prof.var)
tf_data = prof.data[select_tfs.(prof.var.index), :]
tf_s = prof.layers[:Ms][select_tfs.(prof.var.index), :]

tf_data = tf_data[.!(ismissing.(tf_vars.fit_likelihood)), :]
tf_s = tf_s[.!(ismissing.(tf_vars.fit_likelihood)), :]
filter!(:fit_likelihood => x -> !ismissing(x), tf_vars)

## Select target gene

# gene_name = "Rps3"
# i = collect(1:nrow(vars))[vars.index .== gene_name][1]

## Analysis of $s_{tf}$ and $u_{targ}$

# j = 8
# tf_name = tf_vars[j, :index]
# gene_name = vars[i, :index]
# df = DataFrame(X=tf_s[j, :], Y=u[i, :], Cell=prof.obs.clusters, time=prof.obs.latent_time)
# df.logX = log1p.(df.X)
# df.logY = log1p.(df.Y)
# train = hcat(df.logX, df.logY)

# l1 = layer(df, x=:X, y=:Y, color=:Cell, Geom.point)
# p1 = plot(l1,
# 	 Guide.title("Relationship of s_tf and u_targ"),
# 	 Guide.xlabel("Spliced RNA of TF gene, $(tf_name)"),
# 	 Guide.ylabel("Unspliced RNA of target gene, $(gene_name)"),
# )

# l2 = layer(df, x=:logX, y=:logY, color=:Cell, Geom.point)
# p2 = plot(l2,
# 	 Guide.title("Relationship of s_tf and u_targ"),
# 	 Guide.xlabel("log2 spliced RNA of TF gene, $(tf_name)"),
# 	 Guide.ylabel("log2 unspliced RNA of target gene, $(gene_name)"),
# )

# k = 4
# gmm = GaussianMixtures.GMM(k, train, kind=:full)
# posterior, ll = GaussianMixtures.gmmposterior(gmm, train)
# df.cluster = string.(assign_clusters(posterior))
# model = MixtureModel(gmm)
# mix_logpdf(x,y) = logpdf(model, [x,y])

# D = pairwise(Euclidean(), train, dims=1)
# ms = mean(silhouettes(clusters, D))

# xmax = ceil(maximum(df.logX))
# xmin = floor(minimum(df.logX))
# ymax = ceil(maximum(df.logY))
# ymin = floor(minimum(df.logY))

# l3 = layer(df, x=:logX, y=:logY, Geom.point)
# l4 = layer(z=mix_logpdf, xmin=[xmin], xmax=[xmax], ymin=[ymin], ymax=[ymax], Geom.contour)
# coord = Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
# p3 = plot(l3, l4, coord,
#     Guide.xlabel("log spliced RNA of TF gene, $(tf_name)"),
#     Guide.ylabel("log unspliced RNA of target gene, $(gene_name)"),
# )
# p3 |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log likelihood plot.svg"), 10inch, 6inch)


# cluster plot
# l5 = layer(df, x=:logX, y=:logY, color=:cluster, Geom.point,
#     Geom.ellipse(levels=[0.95, 0.99])
# )
# p4 = plot(l5, coord,
#     Guide.xlabel("log spliced RNA of TF gene, $(tf_name)"),
#     Guide.ylabel("log unspliced RNA of target gene, $(gene_name)"),
# )
# p4 |> SVG(joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log cluster plot.svg"), 10inch, 6inch)


# ----------------------------------------------------

function log_likelihood_plot(df, tf_name, gene_name, mix_logpdf;
                             xmax=ceil(maximum(df.logX)), xmin=floor(minimum(df.logX)),
                             ymax = ceil(maximum(df.logY)), ymin = floor(minimum(df.logY)))
    l1 = layer(df, x=:logX, y=:logY, Geom.point)
    l2 = layer(z=mix_logpdf, xmin=[xmin], xmax=[xmax], ymin=[ymin], ymax=[ymax], Geom.contour)
    coord = Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    xlabel = "log spliced RNA of TF gene, $(tf_name)"
    ylabel = "log unspliced RNA of target gene, $(gene_name)"
    filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log likelihood plot.svg")
    plot(l1, l2, coord, Guide.xlabel(xlabel), Guide.ylabel(ylabel)) |> SVG(filepath, 10inch, 6inch)
end

function cluster_plot(df, tf_name, gene_name; xmax=ceil(maximum(df.logX)), xmin=floor(minimum(df.logX)),
    ymax = ceil(maximum(df.logY)), ymin = floor(minimum(df.logY)))
    l = layer(df, x=:logX, y=:logY, color=:clusters, Geom.point)
    coord = Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    xlabel = "log spliced RNA of TF gene, $(tf_name)"
    ylabel = "log unspliced RNA of target gene, $(gene_name)"
    filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log cluster plot.svg")
    plot(l, coord, Guide.xlabel(xlabel), Guide.ylabel(ylabel)) |> SVG(filepath, 10inch, 6inch)
end

function add_logger(func)
    io = open("model-selection.log", "w+")
    logger = SimpleLogger(io)
    with_logger(func, logger)
    close(io)
end

function training_process(k_range, λ, tf_vars, vars; logger=true)
    total_results = []
    splock = Threads.SpinLock()
    process = Threads.@threads for i = 1:nrow(vars)
        for j = 1:nrow(tf_vars)
            tf_name = tf_vars[j, :index]
            gene_name = vars[i, :index]

            df = DataFrame(logX = log1p.(tf_s[j, :]), logY = log1p.(u[i, :]))
            data = hcat(df.logX, df.logY)

            results = grid_search(GMR, data, k_range, λ=λ, verbosity=2)
            best_res = best_result(results; criterion=aic)
            if haskey(best_res, :model)
                best_k = best_res[:k]
                model = best_res[:model]
                scores = best_res[:score]
                clusters = assign_clusters(model, data)
                mix_logpdf(x,y) = logpdf(model.dist, [x,y])

                lock(splock) do
                    @info "(i=$i, j=$j) $tf_name - $gene_name: best k = $best_k"
                    r = (tf_name=tf_name, gene_name=gene_name,
                        best_k=best_k, scores=scores,
                        model=model, clusters=clusters)
                    push!(total_results, r)
                end

                if best_k != 1
                    df.clusters = string.(clusters)
                    log_likelihood_plot(df, tf_name, gene_name, mix_logpdf)
                    cluster_plot(df, tf_name, gene_name)
                end
            end
        end
    end

    logger && add_logger(process)

    total_results
end


k_range = 1:5
λ = 3e-3
# total_results = training_process(k_range, λ, tf_vars, vars)
total_results = []
splock = Threads.SpinLock()
Threads.@threads for i = 1:nrow(vars)
    for j = 1:nrow(tf_vars)
        tf_name = tf_vars[j, :index]
        gene_name = vars[i, :index]

        df = DataFrame(logX = log1p.(tf_s[j, :]), logY = log1p.(u[i, :]))
        data = hcat(df.logX, df.logY)

        results = grid_search(GMR, data, k_range, λ=λ, verbosity=2)
        best_res = best_result(results; criterion=aic)
        if haskey(best_res, :model)
            best_k = best_res[:k]
            model = best_res[:model]
            scores = best_res[:score]
            clusters = assign_clusters(model, data)
            mix_logpdf(x,y) = logpdf(model.dist, [x,y])

            lock(splock) do
                @info "(i=$i, j=$j) $tf_name - $gene_name: best k = $best_k"
                r = (tf_name=tf_name, gene_name=gene_name,
                    best_k=best_k, scores=scores,
                    model=model, clusters=clusters)
                push!(total_results, r)
            end

            if best_k != 1
                df.clusters = string.(clusters)
                log_likelihood_plot(df, tf_name, gene_name, mix_logpdf)
                cluster_plot(df, tf_name, gene_name)
            end
        end
    end
end

# report

report = DataFrame()
report.tf_name = map(x -> x[:tf_name], total_results)
report.gene_name = map(x -> x[:gene_name], total_results)
report.best_k = map(x -> x[:best_k], total_results)
report.scores = map(x -> x[:scores], total_results)
report = report[report.best_k .!= 1, :]
sort!(report, :scores)

save(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results", total_results)
