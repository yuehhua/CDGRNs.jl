using CDGRN
using DataFrames
using FileIO
using JLD2
using SnowyOwl
using Distances
using Clustering
using StatsBase, Statistics
using HypothesisTests
using Gadfly
# import Cairo

function load_profile(dir::String)
    prof = load_data(dir)
    add_unspliced_data!(prof, dir)
    add_velocity!(prof, dir)
    add_moments!(prof, dir)
    return prof
end

function regulation_correlation(filename)
    total_results = load(filename, "total_results")
    nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)
    cor_pairs = corr_table(nonsingle_pairs)
    cor_pairs.is_tf = cor_pairs.tf .== cor_pairs.target
    return cor_pairs, nonsingle_pairs
end

function remove_spurious_pairs(cor_pairs, nonsingle_pairs)
    # map to database
    database = load_CHEA(joinpath(CDGRN.PROJECT_PATH, "CHEA"))
    pairset = CDGRN.make_pairset(database)
    cor_pairs.is_regulation = CDGRN.query_pairset(cor_pairs, pairset)
    true_regulations = cor_pairs[cor_pairs.is_regulation .& .!cor_pairs.is_tf, :]
    true_reg_pairs = filter(x -> (uppercase(x[:tf_name]), uppercase(x[:gene_name])) in pairset, nonsingle_pairs)
    return true_regulations, true_reg_pairs
end

function build_tree(prof, true_reg_pairs; linkage=:ward, save=nothing)
    features = DataFrame(cell=prof.obs.clusters, time=prof.obs.latent_time)
    for res in true_reg_pairs
        colname = res[:tf_name] * "_" * res[:gene_name]
        features[!, colname] = res[:clusters]
    end
    
    data = Array(features[:, 3:end])
    D = pairwise(Hamming(), data, dims=1)
    tree = hclust(D, linkage=linkage, branchorder=:optimal)
    isnothing(save) || CDGRN.clustermap(D, features.cell, filename=save)
    return tree, features
end

function extract_context!(cell_clusters, tree, k)
    cell_clusters[:, Symbol("k$k")] = cutree(tree; k=k)
    return cell_clusters
end

function context_correlation(tfs, prof, true_regulations, context, k)
    pairs = unique(zip(true_regulations.tf, true_regulations.target))
    context_cor = CDGRN.cor(tfs, prof, pairs, context)
    context_pairs = collect(zip(context_cor.tf, context_cor.target))
    context_cor[!, :context], selected_context = CDGRN.max_cor(context_cor, k)
    return context_cor, selected_context, context_pairs
end

function test_pmf(ρ1, ρ2, condition1, condition2; α=0.7, bincount=30, plot_dir=nothing, title="")
    test_result = MannWhitneyUTest(abs.(ρ1), abs.(ρ2))

    if !isnothing(plot_dir)
        ρ = vcat(ρ1, ρ2)
        condition = vcat(repeat([condition1], length(ρ1)), repeat([condition2], length(ρ2)))
        plot_df = DataFrame(ρ=ρ, condition=condition)
        plot_df[!,:alpha] .= α

        Gadfly.with_theme(style(grid_line_style=:solid)) do
            p = plot(
                plot_df, x=:ρ, color=:condition, alpha=:alpha,
                Geom.histogram(position=:identity, bincount=bincount),
                Guide.xlabel("TF-target pair correlation")
            )
            draw(SVG(joinpath(plot_dir, "histogram_$(title).svg"), 6inch, 4inch), p)
            draw(PNG(joinpath(plot_dir, "histogram_$(title).png"), 6inch, 4inch), p)
        end
    end
    return test_result
end

function test_cdf(ρ1, ρ2, condition1, condition2; step=0.1, plot_dir=nothing, title="")
    test_result = ApproximateTwoSampleKSTest(ρ1, ρ2)

    if !isnothing(plot_dir)
        r = -1.0:step:1.0
        cntx_cdf = DataFrame(x=r, y=ecdf(ρ1)(r), condition=condition1)
        global_cdf = DataFrame(x=r, y=ecdf(ρ2)(r), condition=condition2)
        cdf_df = vcat(cntx_cdf, global_cdf)

        Gadfly.with_theme(style(grid_line_style=:solid)) do
            p = plot(
                cdf_df, x=:x, y=:y, color=:condition,
                Geom.step, Guide.xlabel("TF-target pair correlation"),
                Guide.yticks(ticks=[0., 0.25, 0.5, 0.75, 1.])
            )
            draw(SVG(joinpath(plot_dir, "cdf_$(title).svg"), 6inch, 4inch), p)
            draw(PNG(joinpath(plot_dir, "cdf_$(title).png"), 6inch, 4inch), p)
        end
    end
    return test_result
end

## Load data

dir = joinpath(CDGRN.PROJECT_PATH, "results", "pancreas")
prof = load_profile(dir)
tfs = copy(prof)

CDGRN.filter_genes!(prof)
vars = prof.var
u = prof.layers[:Mu]

tf_set = CDGRN.load_tfs(joinpath(dir, "tf_set.jld2"))
CDGRN.filter_tfs!(tfs, tf_set)
tf_vars = tfs.var
tf_s = tfs.layers[:Ms]


filename = joinpath(dir, "GMM-model-selection-result.jld2")
cor_pairs, nonsingle_pairs = regulation_correlation(filename)
true_regulations, true_reg_pairs = remove_spurious_pairs(cor_pairs, nonsingle_pairs)

k = 5
tree, cell_clusters = build_tree(prof, true_reg_pairs, save="clustermap_pancreatic")
extract_context!(cell_clusters, tree, k)

context_cor, selected_context, context_pairs = context_correlation(tfs, prof, true_regulations, cell_clusters.k5, k)
global_ρs = CDGRN.global_cor(tfs, prof, context_pairs, map(x -> x.size, context_cor.context))

# Compare context and global

test_df = innerjoin(context_cor, global_ρs, on=[:tf, :target], makeunique=true)
context_ρ = map(x -> x.ρ, test_df.context)
global_ρ = map(x -> x.ρ, test_df.global)
test_result = test_pmf(context_ρ, global_ρ, "context", "global"; plot_dir="pics/tf-gene gmm model/CDGRN", title="context-global")
test_result = test_cdf(context_ρ, global_ρ, "context", "global"; plot_dir="pics/tf-gene gmm model/CDGRN", title="context-global")


# Compare spliced only

spliced_ρs = CDGRN.spliced_cor(tfs, prof, context_pairs, cell_clusters.k5, selected_context)
test_df2 = innerjoin(context_cor, spliced_ρs, on=[:tf, :target], makeunique=true)
test_df2.spliced_ρ = map(x -> x.ρ, test_df2.spliced)
test_df2.context_ρ = map(x -> x.ρ, test_df2.context)
test_df2 = test_df2[test_df2.tf .!= test_df2.target, :]
test_df2 = test_df2[.!isnan.(test_df2.spliced_ρ), :]
test_result = test_pmf(test_df2.context_ρ, test_df2.spliced_ρ, "unspliced+spliced", "spliced"; plot_dir="pics/tf-gene gmm model/CDGRN", title="unspliced-spliced")
test_result = test_cdf(test_df2.context_ρ, test_df2.spliced_ρ, "unspliced+spliced", "spliced"; plot_dir="pics/tf-gene gmm model/CDGRN", title="unspliced-spliced")


# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)
trainX = Array(df[:, 3:end])'

p = CDGRN.plot_3d_pca(trainX, cell_clusters.k5)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "cell clusters k5.html")
savefig(p, filepath)

# Investigate cell clusters

countmap(cell_clusters[cell_clusters.k3 .== 3, :cell])


# Train on Ngn3 high EP

cortable = Dict()
for c in [1, 3, 4, 5]
    context = cell_clusters.k5 .== c
    cdgrn = train(ContextDependentGRN, tfs, prof, true_regulations, context)
    evaluate!(cdgrn)

    context_cor = CDGRN.cor(cdgrn, tfs, prof, context)
    context_cor = context_cor[.!isnan.(context_cor.ρ), :]
    context_cor.reg_type = context_cor.ρ .> 0
    context_cor.reg_stng = abs.(context_cor.ρ)
    CSV.write(joinpath(dir, "cdgrn_k5", "cdgrn_k5-$(c).csv"), context_cor)

    ## Effective gene set
    correlated = context_cor[context_cor.reg_stng .≥ 0.3, :]
    gene_set = DataFrame(gene=unique!(vcat(correlated.tf, correlated.target)))
    CSV.write(joinpath(dir, "cdgrn_k5", "cdgrn_k5-$(c)_gene_set.csv"), gene_set)

    ## add gene expression
    context_cor.tf_s = [mean(vec(get_gene_expr(tfs, String(g), :Ms))[context]) for g in context_cor.tf]
    context_cor.target_u = [mean(vec(get_gene_expr(prof, String(g), :Mu))[context]) for g in context_cor.target]
    cortable[c] = context_cor
end

## join tables
joined = outerjoin(cortable[1], cortable[3], on=[:tf, :target], renamecols="_cntx1"=>"_cntx3")
joined = outerjoin(joined, cortable[4], on=[:tf, :target], renamecols=""=>"_cntx4")
joined = outerjoin(joined, cortable[5], on=[:tf, :target], renamecols=""=>"_cntx5")
replace!(joined.ρ_cntx1, missing=>0)
replace!(joined.reg_stng_cntx1, missing=>0)
replace!(joined.tf_s_cntx1, missing=>0)
replace!(joined.target_u_cntx1, missing=>0)
replace!(joined.ρ_cntx3, missing=>0)
replace!(joined.reg_stng_cntx3, missing=>0)
replace!(joined.tf_s_cntx3, missing=>0)
replace!(joined.target_u_cntx3, missing=>0)
replace!(joined.ρ_cntx4, missing=>0)
replace!(joined.reg_stng_cntx4, missing=>0)
replace!(joined.tf_s_cntx4, missing=>0)
replace!(joined.target_u_cntx4, missing=>0)
replace!(joined.ρ_cntx5, missing=>0)
replace!(joined.reg_stng_cntx5, missing=>0)
replace!(joined.tf_s_cntx5, missing=>0)
replace!(joined.target_u_cntx5, missing=>0)
CSV.write(joinpath(dir, "cdgrn_k5", "cdgrn_k5_all.csv"), joined)

l3 = Gadfly.layer(context_cor3, x=:ρ, color=[colorant"blue"], Geom.histogram(position=:identity, density=true))
l5 = Gadfly.layer(context_cor5, x=:ρ, color=[colorant"orange"], Geom.histogram(position=:identity, density=true))
l4 = Gadfly.layer(context_cor4, x=:ρ, color=[colorant"red"], Geom.histogram(position=:identity, density=true))
l1 = Gadfly.layer(context_cor1, x=:ρ, color=[colorant"black"], Geom.histogram(position=:identity, density=true))
Gadfly.plot(l3, l5, l4, l1)

sg = CDGRN.build_graph(c)
CDGRN.network_entropy(sg)

# Visualize PCA

df = get_regulation_expr(prof, tfs, true_regulations)
df.context = cell_clusters.k5
trainX = Array(df[:, 3:end])'

p = CDGRN.plot_3d_pca(trainX, df.cell, context)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "cell clusters k3-3.html")
savefig(p, filepath)


# PCA biplot
loading = DataFrame(genes=names(df)[3:end], PC1=model.proj[:,1], PC2=model.proj[:,2], x0=0., y0=0.)
p = Gadfly.plot(loading, x=:x0, y=:y0, xend=:PC1, yend=:PC2, label=:genes,
                Geom.label, Geom.segment(arrow=true))



# Visualize regulation



# k5: 3

target = "Ccne2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "Atad2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "E2f1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-3 regulation $target (global).html")
savefig(p, filepath)


# k5: 2

target = "E2f1"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-2 regulation $target.html")
savefig(p, filepath)


# k5: 4

target = "Pax6"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Pax6", "Pdx1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)

target = "Naaladl2"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Elf5", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Elf5", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Cpe"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Nr3c1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Nr3c1", target)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)


target = "Vdr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target.html")
savefig(p, filepath)
# spliced
p = plot_regulations(data, "Vdr", "E2f1", target, spliced=true)
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (spliced).html")
savefig(p, filepath)
# global
data = make_vis_data(target, tfs, prof, true_regulations)
p = plot_regulations(data, "Vdr", "E2f1", target, model=cdgrn.models[Symbol(target)])
filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "CDGRN", "k5-4 regulation $target (global).html")
savefig(p, filepath)

target = "Rimbp2"
target = "Ghr"
data = make_vis_data(target, tfs, prof, true_regulations, cluster=context)
p = plot_regulations(data, "Elf5", "Pax6", target, model=cdgrn.models[Symbol(target)])
