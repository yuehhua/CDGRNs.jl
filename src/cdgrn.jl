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

function train(::Type{CDGRN}, tfs::Profile, prof::Profile, pairs, context; target::Symbol=:target)
    targets = unique(pairs[!, target])
    models = []
    tf_sets = []
    for targ in targets
        tf_set = unique(pairs.tf[pairs.target .== targ])
        y = vec(get_gene_expr(prof, targ, :Mu))[context]
        X = hcat(map(tf -> vec(get_gene_expr(tfs, tf, :Ms)), tf_set)...)[context, :]
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

function partial_corr(cdgrn::CDGRN, tfs::Profile, prof::Profile, context)
    df = DataFrame(tf=Symbol[], target=Symbol[], ρ=Float64[])
    for targ in cdgrn.scores.target
        tf_set = cdgrn.tfs[targ]
        y = vec(get_gene_expr(prof, string(targ), :Mu))[context]
        xs = map(tf -> vec(get_gene_expr(tfs, string(tf), :Ms)[context]), tf_set)
        ρs = partial_corr(y, xs...)
        for (tf, ρ) in zip(tf_set, ρs)
            push!(df, (tf, targ, ρ))
        end
    end
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

# function plot(grn::CDGRN)
    
# end
