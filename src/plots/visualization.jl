function pca(pc1::AbstractVector, pc2::AbstractVector, pc3::AbstractVector, labels::AbstractVector;
             context=nothing)
    p = if isnothing(context)
        SnowyOwl.Plots.pca(pc1, pc2, pc3; group=labels)
    elseif context isa AbstractVector
        α = map(x -> x ? 1. : 0.2, context)
        colors = [RGB(x...) for x in generate_seaborn_colors(labels)]
        SnowyOwl.Plots.pca(pc1, pc2, pc3; group=labels, markercolor=colors, markeralpha=α)
    else
        throw(ArgumentError("context only accepts vector."))
    end

    return p
end

function pca(pc1::AbstractVector, pc2::AbstractVector, labels::AbstractVector;
             context=nothing, figsize=(1000, 600), dpi=300)
    p = if isnothing(context)
        SnowyOwl.Plots.pca(pc1, pc2; group=labels, figsize=figsize, dpi=dpi)
    elseif context isa AbstractVector
        α = map(x -> x ? 1. : 0.2, context)
        colors = [RGB(x...) for x in generate_seaborn_colors(labels)]
        SnowyOwl.Plots.pca(pc1, pc2; group=labels, markercolor=colors, markeralpha=α,
            figsize=figsize, dpi=dpi)
    else
        throw(ArgumentError("context only accepts vector."))
    end

    return p
end

umap(umap1::AbstractVector, umap2::AbstractVector; figsize=(1000, 600), dpi=300) =
    SnowyOwl.Plots.umap(umap1, umap2; group=:cell, color_palette=:glasbey_hv_n256, figsize=figsize, dpi=dpi)
