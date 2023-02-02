struct ContextDependentGRN{T,S}
    models::Dict{Symbol,T}
    tfs::Dict{Symbol,S}
    scores::DataFrame
end

function train(::Type{ContextDependentGRN}, tfs::OmicsProfile, prof::OmicsProfile, pairs::DataFrame, context; target::Symbol=:target, unspliced=true)
    M = unspliced ? prof.Mu : prof.Ms
    targets = unique(pairs[!, target])
    models = []
    tf_sets = []
    for targ in targets
        tf_set = unique(pairs.tf[pairs.target .== targ])
        y = vec(M[targ, :])[context]
        X = hcat(map(tf -> vec(tfs.Ms[tf, :]), tf_set)...)[context, :]
        push!(models, fit(LinearModel, X, y))
        push!(tf_sets, tf_set)
    end
    res1 = Dict(map(x -> Symbol(x[1])=>x[2], zip(targets, models)))
    res2 = Dict(map(x -> Symbol(x[1])=>Symbol.(x[2]), zip(targets, tf_sets)))
    scores = DataFrame(target=Symbol[], adjR2=Float64[])
    return ContextDependentGRN(res1, res2, scores)
end

function train(::Type{ContextDependentGRN}, tfs::OmicsProfile, prof::OmicsProfile, pairs::DataFrame; target::Symbol=:target)
    targets = unique(pairs[!, target])
    models = []
    tf_sets = []
    for targ in targets
        tf_set = unique(pairs.tf[pairs.target .== targ])
        y = collect(vec(prof.Mu[targ, :]))
        X = hcat(map(tf -> vec(tfs.Ms[tf, :]), tf_set)...)
        push!(models, fit(LinearModel, X, y))
        push!(tf_sets, tf_set)
    end
    res1 = Dict(map(x -> Symbol(x[1])=>x[2], zip(targets, models)))
    res2 = Dict(map(x -> Symbol(x[1])=>Symbol.(x[2]), zip(targets, tf_sets)))
    scores = DataFrame(target=Symbol[], adjR2=Float64[])
    return ContextDependentGRN(res1, res2, scores)
end

function evaluate!(grn::ContextDependentGRN)
    for k in keys(grn.models)
        adjR2 = adjr²(grn.models[k])
        if 0. ≤ adjR2 ≤ 1.0
            push!(grn.scores, (k, adjR2))
        end
    end
    sort!(grn.scores, :adjR2, rev=true)
    return grn
end

function cor(cdgrn::ContextDependentGRN, tfs::OmicsProfile, prof::OmicsProfile, context::BitVector; unspliced=true)
    M = unspliced ? prof.Mu : prof.Ms
    cond = unspliced ? :context : :spliced
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], condition=Symbol[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = vec(M[string(targ), :])[context]
        xs = map(tf -> vec(tfs.Ms[string(tf), :][context]), tf_set)
        ρs = cor(y, xs...)
        for (tf, ρ) in zip(tf_set, ρs)
            if unspliced || tf != targ
                push!(df, (tf, targ, ρ, cond))
            end
        end
    end
    return df
end

function cor(cdgrn::ContextDependentGRN, tfs::OmicsProfile, prof::OmicsProfile; nsample=SnowyOwl.nobs(prof))
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], condition=Symbol[])
    idx = (nsample == SnowyOwl.nobs(prof)) ? Colon() : StatsBase.sample(collect(1:SnowyOwl.nobs(prof)), nsample)
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = collect(vec(prof.Mu[string(targ), :])[idx])
        xs = map(tf -> collect(vec(tfs.Ms[string(tf), :])[idx]), tf_set)
        ρs = cor(y, xs...)
        for (tf, ρ) in zip(tf_set, ρs)
            push!(df, (tf, targ, ρ, :global))
        end
    end
    return df
end

function partial_corr(cdgrn::ContextDependentGRN, tfs::OmicsProfile, prof::OmicsProfile, context; unspliced=true)
    M = unspliced ? prof.Mu : prof.Ms
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], dof=Int[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = vec(M[string(targ), :])[context]
        xs = map(tf -> vec(tfs.Ms[string(tf), :][context]), tf_set)
        ρs = partial_corr(y, xs...)
        dof = length(tf_set) - 1
        for (tf, ρ) in zip(tf_set, ρs)
            if unspliced || tf != targ
                push!(df, (tf, targ, ρ, dof))
            end
        end
    end
    df.dof .= length(df.dof) - 3 .- df.dof
    return df
end

function partial_corr(cdgrn::ContextDependentGRN, tfs::OmicsProfile, prof::OmicsProfile)
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], dof=Int[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = collect(vec(prof.Mu[string(targ), :]))
        xs = map(tf -> collect(vec(tfs.Ms[string(tf), :])), tf_set)
        ρs = partial_corr(y, xs...)
        dof = length(tf_set) - 1
        for (tf, ρ) in zip(tf_set, ρs)
            push!(df, (tf, targ, ρ, dof))
        end
    end
    df.dof .= length(df.dof) - 3 .- df.dof
    return df
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
        context_cor.tf_s = [mean(vec(tfs.Ms[String(g), :])[cntx]) for g in context_cor.tf]
        context_cor.target_u = [mean(vec(prof.Mu[String(g), :])[cntx]) for g in context_cor.target]
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
