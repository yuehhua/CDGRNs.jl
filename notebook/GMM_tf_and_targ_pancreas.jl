using Logging

using CDGRN
using SnowyOwl
using DataFrames
using JLD2
using FileIO
using Gadfly
using Statistics
Gadfly.set_default_plot_size(8inch, 6inch)

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "pancreas")
prof = load_profile(dir)
tf_set = CDGRN.load_tfs(joinpath(dir, "tf_set.jld2"))
tfs = select_genes!(copy(prof), tf_set)

select_high_likelihood!(prof)
vars = prof.var
u = prof.layers[:Mu]

select_high_likelihood!(tfs, min_likelihood=-Inf)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]

function log_likelihood_plot(df, tf_name, gene_name, mix_logpdf;
                             xmax=ceil(maximum(df.logX)), xmin=floor(minimum(df.logX)),
                             ymax = ceil(maximum(df.logY)), ymin = floor(minimum(df.logY)))
    l1 = layer(df, x=:logX, y=:logY, Geom.point)
    l2 = layer(z=mix_logpdf, xmin=[xmin], xmax=[xmax], ymin=[ymin], ymax=[ymax], Geom.contour)
    coord = Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    xlabel = "log spliced RNA of TF gene, $(tf_name)"
    ylabel = "log unspliced RNA of target gene, $(gene_name)"
    filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log likelihood plot.svg")
    plot(l1, l2, coord, Guide.xlabel(xlabel), Guide.ylabel(ylabel)) |> SVG(filepath, 10inch, 6inch)
end

function cluster_plot(df, tf_name, gene_name; xmax=ceil(maximum(df.logX)), xmin=floor(minimum(df.logX)),
    ymax = ceil(maximum(df.logY)), ymin = floor(minimum(df.logY)))
    l = layer(df, x=:logX, y=:logY, color=:clusters, Geom.point)
    coord = Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    xlabel = "log spliced RNA of TF gene, $(tf_name)"
    ylabel = "log unspliced RNA of target gene, $(gene_name)"
    filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "$(tf_name)-$(gene_name) log cluster plot.svg")
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

failed = Channel(Inf)

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

        try
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
        catch e
            put!(failed, (tf_name, gene_name))
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
