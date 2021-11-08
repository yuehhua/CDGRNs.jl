using CDGRN
using JLD2
using Gadfly

@load "results/model-selection-result.jld2" all_pairs

p = plot(all_pairs, x=:ll, color=:best_k, Geom.histogram,
         Guide.xlabel("Log Likelihood"))
p |> SVG(joinpath(CDGRN.PROJECT_PATH, "pics", "model-selection", "log likelihood distribution plot-total.svg"), 10inch, 6inch)


p2 = plot(all_pairs[all_pairs.best_k .> 1, :], x=:ll, color=:best_k, Geom.histogram,
          Guide.xlabel("Log Likelihood"))
p2 |> SVG(joinpath(CDGRN.PROJECT_PATH, "pics", "model-selection", "log likelihood distribution plot-k:2~5.svg"), 10inch, 6inch)


p3 = plot(all_pairs[all_pairs.best_k .> 2, :], x=:ll, color=:best_k, Geom.histogram,
          Guide.xlabel("Log Likelihood"))
p3 |> SVG(joinpath(CDGRN.PROJECT_PATH, "pics", "model-selection", "log likelihood distribution plot-k:3~5.svg"), 10inch, 6inch)

