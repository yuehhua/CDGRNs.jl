using CDGRN
using DataFrames
using FileIO
using JLD2
using Gadfly
using SnowyOwl
using Distances
Gadfly.set_default_plot_size(8inch, 6inch)
import Cairo

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)
add_moments!(prof, dir)


# correlation analysis

total_results = load(joinpath(dir, "GMM-model-selection-result.jld2"), "total_results")
nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)
correlations = map(x -> correlation(x[:model]), nonsingle_pairs)
corr_dist = CDGRN.fisher_transform(vcat(correlations...))

p = plot(x=corr_dist, Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist.svg")
p |> SVG(filepath, 5inch, 3inch)


# map to curated TF-target database

## unit: component
cor_pairs = map(x -> (tf=x[:tf_name], target=x[:gene_name], corr=correlation(x[:model])), nonsingle_pairs)
cor_pairs = flatten(DataFrame(cor_pairs), :corr)
cor_pairs.adjusted_corr = CDGRN.fisher_transform(cor_pairs.corr)

reg_pairs = Set(map(i -> (regulations[i,:tf], regulations[i,:target]), 1:nrow(regulations)))

query_pairs = map(i -> (uppercase(cor_pairs[i,:tf]), uppercase(cor_pairs[i,:target])), 1:nrow(cor_pairs))
cor_pairs.is_regulation = map(x -> x in reg_pairs, query_pairs)  # Components: 2270

unique(cor_pairs[cor_pairs.is_regulation, :], [:tf, :target])  # TF-gene pairs: 830

p = plot(cor_pairs, x=:adjusted_corr, color=:is_regulation,
         Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist with label.svg")
p |> SVG(filepath, 5inch, 3inch)

p = plot(x=cor_pairs[cor_pairs.is_regulation, :adjusted_corr],
         Geom.histogram, Guide.xlabel("TF-gene pair component correlation (Fisher's z-transformed)"))
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "correlation hist (true regulation).svg")
p |> SVG(filepath, 5inch, 3inch)

## unit: pair
cor_pairs = map(x -> (tf=x[:tf_name], target=x[:gene_name], corr=correlation(x[:model])), nonsingle_pairs)
cor_pairs = DataFrame(cor_pairs)
# The maximum absolute (component) correlation as pair correlation
cor_pairs.corr = map(x -> x[argmax(abs.(x))], cor_pairs.corr)
cor_pairs.adjusted_corr = CDGRN.fisher_transform(cor_pairs.corr)

using HypothesisTests

query_pairs = map(i -> (uppercase(cor_pairs[i,:tf]), uppercase(cor_pairs[i,:target])), 1:nrow(cor_pairs))
cor_pairs.is_regulation = map(x -> x in reg_pairs, query_pairs)

tests = []
for z in 2:150
    xs = cor_pairs.adjusted_corr .> z
    ys = cor_pairs.is_regulation
    a = count(xs .& ys)
    b = count(.!(xs) .& ys)
    c = count(xs .& .!(ys))
    d = count(.!(xs) .& .!(ys))
    test = HypothesisTests.FisherExactTest(a, b, c, d)
    push!(tests, test)
end
pvalues = pvalue.(tests)
p = plot(x=2:150, y=pvalues, Geom.line,
         Guide.xlabel("Fisher transformed z"), Guide.ylabel("p value"))
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "p-value.svg")
p |> SVG(filepath, 5inch, 3inch)

scores = -log10.(pvalues)
p = plot(x=2:150, y=scores, Geom.line,
         Guide.xlabel("Fisher transformed z"), Guide.ylabel("Score(- log 10 (p value))"))
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "correlations", "score.svg")
p |> SVG(filepath, 5inch, 3inch)
