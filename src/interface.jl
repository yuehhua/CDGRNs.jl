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
    database = load_CHEA(joinpath(PROJECT_PATH, "CHEA"))
    pairset = make_pairset(database)
    cor_pairs.is_regulation = query_pairset(cor_pairs, pairset)
    true_regulations = cor_pairs[cor_pairs.is_regulation .& .!cor_pairs.is_tf, :]
    true_reg_pairs = filter(x -> (uppercase(x[:tf_name]), uppercase(x[:gene_name])) in pairset, nonsingle_pairs)
    return true_regulations, true_reg_pairs
end

function build_tree(prof, true_reg_pairs; linkage=:ward, celltype=:clusters,
                    time=:latent_time, col_palette=:default, save=nothing)
    features = DataFrame()
    isnothing(celltype) || (features[!, :cell] = prof.obs[!, celltype])
    isnothing(time) || (features[!, :time] = prof.obs[!, time])
    for res in true_reg_pairs
        colname = res[:tf_name] * "_" * res[:gene_name]
        features[!, colname] = res[:clusters]
    end
    
    data = Array(features[:, 3:end])
    D = pairwise(Hamming(), data, dims=1)
    tree = hclust(D, linkage=linkage, branchorder=:optimal)
    isnothing(save) || clustermap(D, features.cell, save, col_palette=col_palette)
    return tree, features
end

function extract_context!(cell_clusters, tree, k)
    cell_clusters[:, Symbol("k$k")] = cutree(tree; k=k)
    return cell_clusters
end

function context_correlation(tfs, prof, true_regulations, context, k)
    pairs = unique(zip(true_regulations.tf, true_regulations.target))
    context_cor = cor(tfs, prof, pairs, context)
    context_pairs = collect(zip(context_cor.tf, context_cor.target))
    context_cor[!, :context], selected_context = max_cor(context_cor, k)
    return context_cor, selected_context, context_pairs
end

function test_pmf(ρ1, ρ2, condition1, condition2; α=0.7, bincount=30,
        plot_dir=nothing, title="", save_png=true, save_svg=true, figsize=(800, 600))
    test_result = MannWhitneyUTest(abs.(ρ1), abs.(ρ2))

    if !isnothing(plot_dir)
        ρ = vcat(ρ1, ρ2)
        condition = vcat(repeat([condition1], length(ρ1)), repeat([condition2], length(ρ2)))
        plot_df = DataFrame(ρ=ρ, condition=condition)
        bins = range(minimum(plot_df[!, :ρ]), stop=maximum(plot_df[!, :ρ]), length=bincount)

        p = @df plot_df histogram(:ρ, group=:condition, bins=bins,
            fillalpha=α, bar_edges=false,
            xlabel="TF-target pair correlation", ylabel="Count",
            thickness_scaling=2, widen=false, size=figsize,
        )
        save_svg && savefig(p, joinpath(plot_dir, "histogram_$(title).svg"))
        save_png && savefig(p, joinpath(plot_dir, "histogram_$(title).png"))
    end
    return test_result
end

function test_cdf(ρ1, ρ2, condition1, condition2; step=0.1,
        plot_dir=nothing, title="", save_png=true, save_svg=true, figsize=(800, 600))
    test_result = ApproximateTwoSampleKSTest(ρ1, ρ2)

    if !isnothing(plot_dir)
        r = -1.0:step:1.0
        cntx_cdf = DataFrame(x=r, y=ecdf(ρ1)(r), condition=condition1)
        global_cdf = DataFrame(x=r, y=ecdf(ρ2)(r), condition=condition2)
        cdf_df = vcat(cntx_cdf, global_cdf)

        default(size = (800, 600))
        p = @df cdf_df Plots.plot(:x, :y, group=:condition, linetype=:steppre,
            xlabel="TF-target pair correlation", ylabel="Cumulative count",
            yticks=[0., 0.25, 0.5, 0.75, 1.], legend_position=:bottomright,
            thickness_scaling=2, widen=false, size=figsize,
        )
        save_svg && savefig(p, joinpath(plot_dir, "cdf_$(title).svg"))
        save_png && savefig(p, joinpath(plot_dir, "cdf_$(title).png"))
    end
    return test_result
end

function save_cdgrn(filename, cdgrn, tfs, prof, cntx)
    context_cor = CDGRNs.cor(cdgrn, tfs, prof, cntx)
    context_cor = context_cor[.!isnan.(context_cor.ρ), :]
    context_cor.reg_type = context_cor.ρ .> 0
    context_cor.reg_stng = abs.(context_cor.ρ)
    CSV.write(filename, context_cor)
    return context_cor
end

function save_effective_gene_set(filename, context_cor, reg_strength)
    correlated = context_cor[context_cor.reg_stng .≥ reg_strength, :]
    gene_set = DataFrame(gene=unique!(vcat(correlated.tf, correlated.target)))
    CSV.write(filename, gene_set)
end

function train_cdgrns(tfs, prof, true_regulations, clusters, selected::AbstractVector{T},
                      reg_strength::Real, dir::String, prefix::String;
                      return_model=false) where {T}
    cortable = Dict{T,DataFrame}()
    cdgrns = Dict{T,ContextDependentGRN}()
    for i in selected
        cntx = clusters .== i
        cdgrn = train(ContextDependentGRN, tfs, prof, true_regulations, cntx)
        evaluate!(cdgrn)
    
        filename = joinpath(dir, "$(prefix)-$(i).csv")
        context_cor = save_cdgrn(filename, cdgrn, tfs, prof, cntx)
        
        filename = joinpath(dir, "$(prefix)-$(i)_gene_set.csv")
        save_effective_gene_set(filename, context_cor, reg_strength)
    
        ## add gene expression
        context_cor.tf_s = [mean(vec(get_gene_expr(tfs, String(g), :Ms))[cntx]) for g in context_cor.tf]
        context_cor.target_u = [mean(vec(get_gene_expr(prof, String(g), :Mu))[cntx]) for g in context_cor.target]
        cortable[i] = context_cor
        return_model && (cdgrns[i] = cdgrn)
    end
    return return_model ? (cortable, cdgrns) : cortable
end

function train_cdgrns(tfs, prof, true_regulations, contexts, selected_contexts::AbstractVector{Int},
                      dir::String; reg_strength=0.3, prefix="cdgrn_k",
                      return_model=false)
    if return_model
        cortable, cdgrns = train_cdgrns(tfs, prof, true_regulations, contexts, selected_contexts,
            reg_strength, dir, prefix, return_model=return_model)
    else
        cortable = train_cdgrns(tfs, prof, true_regulations, contexts, selected_contexts,
            reg_strength, dir, prefix, return_model=return_model)
    end

    i, j = selected_contexts[1], selected_contexts[2]
    joined = outerjoin(cortable[i], cortable[j], on=[:tf, :target], renamecols="_cntx$i"=>"_cntx$j")
    for i in selected_contexts[3:end]
        joined = outerjoin(joined, cortable[i], on=[:tf, :target], renamecols=""=>"_cntx$i")
    end

    for i in selected_contexts
        replace!(joined[!, Symbol("ρ_cntx$i")], missing=>0)
        replace!(joined[!, Symbol("reg_stng_cntx$i")], missing=>0)
        replace!(joined[!, Symbol("tf_s_cntx$i")], missing=>0)
        replace!(joined[!, Symbol("target_u_cntx$i")], missing=>0)
    end
    CSV.write(joinpath(dir, "$(prefix)_all.csv"), joined)
    
    return return_model ? (cortable, cdgrns) : cortable
end

function train_cdgrns(tfs, prof, true_regulations, cells, selected_cells::AbstractVector{String},
                      dir::String; reg_strength=0.3, prefix="cdgrn_k")
    return train_cdgrns(tfs, prof, true_regulations, cells, selected_cells,
                        reg_strength, dir, prefix)
end
