using Logging

using CDGRNs
using GaussianMixtureRegressions
using SnowyOwl
using DataFrames
using JLD2
using FileIO
using Statistics

## Load data

dir = joinpath(CDGRNs.PROJECT_PATH, "results", "pancreas")
fig_dir = joinpath(CDGRNs.PROJECT_PATH, "pics", "pancreas")
prof = CDGRNs.Preprocess.load_profile(dir)
tf_set = CDGRNs.Preprocess.load_tfs(joinpath(dir, "tf_set.jld2"))
tfs = CDGRNs.Preprocess.select_genes!(copy(prof), tf_set)

CDGRNs.Preprocess.select_high_likelihood!(prof)
vars = prof.var
u = prof.layers[:Mu]

CDGRNs.Preprocess.select_high_likelihood!(tfs, min_likelihood=-Inf)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]

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

            try
                results = GaussianMixtureRegressions.grid_search(GaussianMixtureRegressions.GMR, data, k_range, λ=λ, verbosity=2)
                best_res = GaussianMixtureRegressions.best_result(results; criterion=GaussianMixtureRegressions.aic)
                if haskey(best_res, :model)
                    best_k = best_res[:k]
                    model = best_res[:model]
                    scores = best_res[:score]
                    clusters = GaussianMixtureRegressions.assign_clusters(model, data)
                    mix_logpdf(x,y) = logpdf(model.dist, [x,y])
    
                    lock(splock) do
                        @info "(i=$i, j=$j) $tf_name - $gene_name: best k = $best_k"
                        r = (tf_name=tf_name, gene_name=gene_name,
                            best_k=best_k, scores=scores,
                            model=model, clusters=clusters)
                        push!(total_results, r)
                    end
    
                    if best_k != 1
                        CDGRNs.Plots.likelihood_landscape(
                            df, tf_name, gene_name, mix_logpdf, savepath=fig_dir
                        )
                        CDGRNs.Plots.cluster_landscape(
                            df, tf_name, gene_name, string.(clusters), savepath=fig_dir
                        )
                    end
                end
            catch e
                put!(failed, (tf_name, gene_name))
            end
        end
    end

    logger && add_logger(process)

    total_results
end

failed = Channel(Inf)

k_range = 1:5
λ = 3e-3
total_results = training_process(k_range, λ, tf_vars, vars)

# report

report = DataFrame()
report.tf_name = map(x -> x[:tf_name], total_results)
report.gene_name = map(x -> x[:gene_name], total_results)
report.best_k = map(x -> x[:best_k], total_results)
report.scores = map(x -> x[:scores], total_results)
report = report[report.best_k .!= 1, :]
sort!(report, :scores)

save(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results", total_results)
