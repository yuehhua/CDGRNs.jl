function context_correlation(tfs, prof, true_regulations, context, k)
    pairs = unique(zip(true_regulations.tf, true_regulations.target))
    context_cor = cor(tfs, prof, pairs, context)
    context_pairs = collect(zip(context_cor.tf, context_cor.target))
    context_cor[!, :context], selected_context = max_cor(context_cor, k)
    return context_cor, selected_context, context_pairs
end

function cor(tfs::OmicsProfile, prof::OmicsProfile, pairs::Vector, clusters::Vector{<:Integer}; unspliced=true)
    M = unspliced ? :Mu : :Ms
    df = DataFrame()
    df[!, :tf] = map(x -> x[1], pairs)
    df[!, :target] = map(x -> x[2], pairs)
    nclst = maximum(clusters)
    not_nan = trues(length(pairs))
    for n in 1:nclst
        context = clusters .== n
        ρs = Float64[]
        for (tf, targ) in pairs
            y = vec(get_gene_expr(prof, targ, M))[context]
            x = vec(get_gene_expr(tfs, tf, :Ms))[context]
            push!(ρs, Statistics.cor(y, x))
        end
        df[!, Symbol("context$(n)")] = map(x -> (ρ=x, size=count(context)), ρs)
        not_nan .&= .!isnan.(ρs)
    end

    return df[not_nan, :]
end

function global_cor(tfs::OmicsProfile, prof::OmicsProfile, pairs::Vector, nsample::Vector{<:Integer})
    df = DataFrame()
    df[!, :tf] = map(x -> x[1], pairs)
    df[!, :target] = map(x -> x[2], pairs)
    res = NamedTuple[]
    for i in 1:length(pairs)
        tf, targ = pairs[i]
        idx = StatsBase.sample(collect(1:SnowyOwl.nobs(prof)), nsample[i])
        y = vec(get_gene_expr(prof, targ, :Mu))[idx]
        x = vec(get_gene_expr(tfs, tf, :Ms))[idx]
        push!(res, (ρ=Statistics.cor(y, x), size=nsample[i]))
    end
    df[!, :global] = res
    return df[map(x -> !isnan(x.ρ), res), :]
end

function spliced_cor(tfs::OmicsProfile, prof::OmicsProfile, pairs::Vector, clusters::Vector{<:Integer}, context::Vector{<:Integer})
    df = DataFrame()
    df[!, :tf] = map(x -> x[1], pairs)
    df[!, :target] = map(x -> x[2], pairs)
    res = NamedTuple[]
    for i in 1:length(pairs)
        tf, targ = pairs[i]
        idx = clusters .== context[i]
        y = vec(get_gene_expr(prof, targ, :Ms))[idx]
        x = vec(get_gene_expr(tfs, tf, :Ms))[idx]
        push!(res, (ρ=Statistics.cor(y, x), size=count(idx)))
    end
    df[!, :spliced] = res
    return df
end

function max_cor(context_cor::DataFrame, nclst::Int)
    context_ρs = []
    selected_context = Int[]
    for i in 1:nrow(context_cor)
        entity = [context_cor[i, Symbol("context$k")] for k in 1:nclst]
        j = argmax(map(x -> abs(x.ρ), entity))
        push!(context_ρs, entity[j])
        push!(selected_context, j)
    end
    return context_ρs, selected_context
end

test_pmf(ρ1, ρ2) = MannWhitneyUTest(abs.(ρ1), abs.(ρ2))

test_cdf(ρ1, ρ2) = ApproximateTwoSampleKSTest(ρ1, ρ2)
