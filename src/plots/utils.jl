function generate_plots_colors(xs::AbstractVector; col_palette=:default)
    keys = sort(unique(xs))
    lut = Dict(zip(keys, palette(col_palette)))
    lut = Dict((k, (v.r, v.g, v.b)) for (k, v) in lut)
    return map(x -> lut[x], xs)
end
