cor(x::AbstractVector, y::AbstractVector, zs...) = [Statistics.cor(x, y), cor(x, zs...)...]
cor(x::AxisArray, y::AxisArray) = Statistics.cor(collect(x), collect(y))

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

cor2dist(ρ::Real) = sqrt(2*(1-ρ))
