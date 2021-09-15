function corr_table(pairs)
    # unit: component
    cor_pairs = map(x -> (tf=x[:tf_name], target=x[:gene_name], best_k=x[:best_k], corr=correlation(x[:model])), pairs)
    df = DataFrames.flatten(DataFrame(cor_pairs), :corr)
    df.adjusted_corr = GRN.fisher_transform(df.corr)
    return df
end

function make_pairset(reg::DataFrame)
    return Set(map((x,y) -> (uppercase(x), uppercase(y)), reg.tf, reg.target))
end

function query_pairset(corr_pairs::DataFrame, reg_pairs::Set)
    query_pairs = map((x,y) -> (uppercase(x), uppercase(y)), corr_pairs.tf, corr_pairs.target)
    return map(x -> x in reg_pairs, query_pairs)
end

struct CDGRN{T,S}
    models::Dict{Symbol,T}
    tfs::Dict{Symbol,S}
    scores::DataFrame
end

function train(::Type{CDGRN}, tfs::Profile, prof::Profile, pairs::DataFrame, context; target::Symbol=:target, unspliced=true)
    M = unspliced ? :Mu : :Ms
    targets = unique(pairs[!, target])
    models = []
    tf_sets = []
    for targ in targets
        tf_set = unique(pairs.tf[pairs.target .== targ])
        y = vec(get_gene_expr(prof, targ, M))[context]
        X = hcat(map(tf -> vec(get_gene_expr(tfs, tf, :Ms)), tf_set)...)[context, :]
        push!(models, fit(LinearModel, X, y))
        push!(tf_sets, tf_set)
    end
    res1 = Dict(map(x -> Symbol(x[1])=>x[2], zip(targets, models)))
    res2 = Dict(map(x -> Symbol(x[1])=>Symbol.(x[2]), zip(targets, tf_sets)))
    scores = DataFrame(target=Symbol[], adjR2=Float64[])
    return CDGRN(res1, res2, scores)
end

function train(::Type{CDGRN}, tfs::Profile, prof::Profile, pairs::DataFrame; target::Symbol=:target)
    targets = unique(pairs[!, target])
    models = []
    tf_sets = []
    for targ in targets
        tf_set = unique(pairs.tf[pairs.target .== targ])
        y = collect(vec(get_gene_expr(prof, targ, :Mu)))
        X = hcat(map(tf -> vec(get_gene_expr(tfs, tf, :Ms)), tf_set)...)
        push!(models, fit(LinearModel, X, y))
        push!(tf_sets, tf_set)
    end
    res1 = Dict(map(x -> Symbol(x[1])=>x[2], zip(targets, models)))
    res2 = Dict(map(x -> Symbol(x[1])=>Symbol.(x[2]), zip(targets, tf_sets)))
    scores = DataFrame(target=Symbol[], adjR2=Float64[])
    return CDGRN(res1, res2, scores)
end

function evaluate!(grn::CDGRN)
    for k in keys(grn.models)
        adjR2 = adjr²(grn.models[k])
        if 0. ≤ adjR2 ≤ 1.0
            push!(grn.scores, (k, adjR2))
        end
    end
    sort!(grn.scores, :adjR2, rev=true)
    return grn
end

function cor(tfs::Profile, prof::Profile, pairs::Vector, clusters::Vector{<:Integer}; unspliced=true)
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
            push!(ρs, cor(y, x))
        end
        df[!, Symbol("context$(n)")] = map(x -> (ρ=x, size=count(context)), ρs)
        not_nan .&= .!isnan.(ρs)
    end

    return df[not_nan, :]
end

function cor(cdgrn::CDGRN, tfs::Profile, prof::Profile, context::BitVector; unspliced=true)
    M = unspliced ? :Mu : :Ms
    cond = unspliced ? :context : :spliced
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], condition=Symbol[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = vec(get_gene_expr(prof, string(targ), M))[context]
        xs = map(tf -> vec(get_gene_expr(tfs, string(tf), :Ms)[context]), tf_set)
        ρs = cor(y, xs...)
        for (tf, ρ) in zip(tf_set, ρs)
            if unspliced || tf != targ
                push!(df, (tf, targ, ρ, cond))
            end
        end
    end
    return df
end

function cor(cdgrn::CDGRN, tfs::Profile, prof::Profile; nsample=SnowyOwl.nobs(prof))
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], condition=Symbol[])
    idx = (nsample == SnowyOwl.nobs(prof)) ? Colon() : StatsBase.sample(collect(1:SnowyOwl.nobs(prof)), nsample)
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = collect(vec(get_gene_expr(prof, string(targ), :Mu))[idx])
        xs = map(tf -> collect(vec(get_gene_expr(tfs, string(tf), :Ms))[idx]), tf_set)
        ρs = cor(y, xs...)
        for (tf, ρ) in zip(tf_set, ρs)
            push!(df, (tf, targ, ρ, :global))
        end
    end
    return df
end

function global_cor(tfs::Profile, prof::Profile, pairs::Vector, nsample::Vector{<:Integer})
    df = DataFrame()
    df[!, :tf] = map(x -> x[1], pairs)
    df[!, :target] = map(x -> x[2], pairs)
    res = NamedTuple[]
    for i in 1:length(pairs)
        tf, targ = pairs[i]
        idx = StatsBase.sample(collect(1:SnowyOwl.nobs(prof)), nsample[i])
        y = vec(get_gene_expr(prof, targ, :Mu))[idx]
        x = vec(get_gene_expr(tfs, tf, :Ms))[idx]
        push!(res, (ρ=cor(y, x), size=nsample[i]))
    end
    df[!, :global] = res
    return df[map(x -> !isnan(x.ρ), res), :]
end

function spliced_cor(tfs::Profile, prof::Profile, pairs::Vector, clusters::Vector{<:Integer}, context::Vector{<:Integer})
    df = DataFrame()
    df[!, :tf] = map(x -> x[1], pairs)
    df[!, :target] = map(x -> x[2], pairs)
    res = NamedTuple[]
    for i in 1:length(pairs)
        tf, targ = pairs[i]
        idx = clusters .== context[i]
        y = vec(get_gene_expr(prof, targ, :Ms))[idx]
        x = vec(get_gene_expr(tfs, tf, :Ms))[idx]
        push!(res, (ρ=cor(y, x), size=count(idx)))
    end
    df[!, :spliced] = res
    return df
end

cor(x::AbstractVector, y::AbstractVector, zs...) = [cor(x, y), cor(x, zs...)...]

function partial_corr(cdgrn::CDGRN, tfs::Profile, prof::Profile, context; unspliced=true)
    M = unspliced ? :Mu : :Ms
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], dof=Int[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = vec(get_gene_expr(prof, string(targ), M))[context]
        xs = map(tf -> vec(get_gene_expr(tfs, string(tf), :Ms)[context]), tf_set)
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

function partial_corr(cdgrn::CDGRN, tfs::Profile, prof::Profile)
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[], dof=Int[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = collect(vec(get_gene_expr(prof, string(targ), :Mu)))
        xs = map(tf -> collect(vec(get_gene_expr(tfs, string(tf), :Ms))), tf_set)
        ρs = partial_corr(y, xs...)
        dof = length(tf_set) - 1
        for (tf, ρ) in zip(tf_set, ρs)
            push!(df, (tf, targ, ρ, dof))
        end
    end
    df.dof .= length(df.dof) - 3 .- df.dof
    return df
end

partial_corr(xs::AbstractVector, ys::AbstractVector) = cor(xs, ys)

function partial_corr(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector)
    ρxy = cor(xs, ys)
    ρxz = cor(xs, zs)
    ρyz = cor(ys, zs)
    ρxy_z = partial_corr(ρxy, ρxz, ρyz)
    ρxz_y = partial_corr(ρxz, ρxy, ρyz)
    return ρxy_z, ρxz_y
end

"""
Return ρxy|z
"""
function partial_corr(ρxy::T, ρxz::T, ρyz::T) where {T<:Real}
    num = ρxy - ρxz * ρyz
    denom = sqrt(1 - ρxz^2) * sqrt(1 - ρyz^2)
    return num / denom
end

function partial_corr(ws::AbstractVector, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector)
    ρwx = cor(ws, xs)
    ρwy = cor(ws, ys)
    ρwz = cor(ws, zs)
    ρxz = cor(xs, zs)
    ρyz = cor(ys, zs)

    ρwx_y, ρwy_x = partial_corr(ws, xs, ys)
    ρwz_x = partial_corr(ρwz, ρwx, ρxz)
    ρwz_y = partial_corr(ρwz, ρwy, ρyz)
    ρxz_y, ρyz_x = partial_corr(zs, xs, ys)

    ρwx_yz = partial_corr(ρwx_y, ρwz_y, ρxz_y)
    ρwy_xz = partial_corr(ρwy_x, ρwz_x, ρyz_x)
    ρwz_xy = partial_corr(ρwz_x, ρwy_x, ρyz_x)
    return ρwx_yz, ρwy_xz, ρwz_xy
end

function partial_corr(vs::AbstractVector, ws::AbstractVector, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector)
    ρvw = cor(vs, ws)
    ρvx = cor(vs, xs)
    ρvy = cor(vs, ys)
    ρvz = cor(vs, zs)
    ρwx = cor(ws, xs)
    ρwy = cor(ws, ys)
    ρwz = cor(ws, zs)
    ρxy = cor(xs, ys)
    ρxz = cor(xs, zs)
    ρyz = cor(ys, zs)

    ρvx_w = partial_corr(ρvx, ρvw, ρwx)
    ρvy_w = partial_corr(ρvy, ρvw, ρwy)
    ρvy_x = partial_corr(ρvy, ρvx, ρxy)
    ρvz_w = partial_corr(ρvz, ρvw, ρwz)
    ρvz_x = partial_corr(ρvz, ρvx, ρxz)
    ρxz_w = partial_corr(ρxz, ρwx, ρwz)
    ρyz_w = partial_corr(ρyz, ρwy, ρwz)
    ρyz_x = partial_corr(ρyz, ρxy, ρxz)

    ρvw_xy, ρvx_wy, ρvy_wx = partial_corr(vs, ws, xs, ys)
    ρvz_wx = partial_corr(ρvz_w, ρvx_w, ρxz_w)
    ρvz_wy = partial_corr(ρvz_w, ρvy_w, ρyz_w)
    ρvz_xy = partial_corr(ρvz_x, ρvy_x, ρyz_x)
    ρwz_xy, ρxz_wy, ρyz_wx = partial_corr(zs, ws, xs, ys)

    ρvw_xyz = partial_corr(ρvw_xy, ρvz_xy, ρwz_xy)
    ρvx_wyz = partial_corr(ρvx_wy, ρvz_wy, ρxz_wy)
    ρvy_wxz = partial_corr(ρvy_wx, ρvz_wx, ρyz_wx)
    ρvz_wxy = partial_corr(ρvz_wx, ρvy_wx, ρyz_wx)
    return ρvw_xyz, ρvx_wyz, ρvy_wxz, ρvz_wxy
end

function to_graph(part_cor::DataFrame; src=:tf, dst=:target,
                  nodes=unique!(vcat(part_cor[:,src], part_cor.target[:,dst])))
    g = MetaDiGraph(length(nodes))
    for row in eachrow(part_cor)
        i = findfirst(nodes .== row[src])
        j = findfirst(nodes .== row[dst])
        add_edge!(g, i, j, :edge_width, row.ρ)
    end
    return g
end

# function network_entropy(g::AbstractGraph)
#     d = degree(g)
#     N = length(d)
#     return sum(log, d) / (N * log(N - 1))
# end
