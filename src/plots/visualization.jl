function pca(pc1::AbstractVector, pc2::AbstractVector, pc3::AbstractVector, labels::AbstractVector;
             context=nothing, savepath::String=@__DIR__, title="", save=[:html])
    p = if isnothing(context)
        Plots.scatter(pc1, pc2, pc3, group=labels, xlabel="PC1", ylabel="PC2", zlabel="PC3",
            markersize=1, markerstrokewidth=0)
    elseif context isa AbstractVector
        α = map(x -> x ? 1. : 0.2, context)
        colors = [RGB(x...) for x in generate_seaborn_colors(labels)]
        Plots.scatter(pc1, pc2, pc3, group=labels, xlabel="PC1", ylabel="PC2", zlabel="PC3",
            markersize=1, markerstrokewidth=0, markercolor=colors, markeralpha=α)
    else
        throw(ArgumentError("context only accepts vector."))
    end
    savefig(p, joinpath(savepath, "pca_$(title)"), save)
    return p
end

function pca(pc1::AbstractVector, pc2::AbstractVector, labels::AbstractVector;
             context=nothing, savepath::String=@__DIR__, title="", save=[:svg, :png],
             figsize=(1000, 600), dpi=300)
    p = if isnothing(context)
        Plots.scatter(pc1, pc2, group=labels, xlabel="PC1", ylabel="PC2",
        markercolor=colors, markeralpha=α,
        legend=:outerright, markersize=2, markerstrokewidth=0,
        thickness_scaling=2, widen=false,
        size=figsize, dpi=dpi)
    elseif context isa AbstractVector
        α = map(x -> x ? 1. : 0.2, context)
        colors = [RGB(x...) for x in generate_seaborn_colors(labels)]
        Plots.scatter(pc1, pc2, group=labels, markercolor=colors, markeralpha=α,
                        xlabel="PC1", ylabel="PC2",
                        legend=:outerright, markersize=2, markerstrokewidth=0,
                        thickness_scaling=2, widen=false,
                        size=figsize, dpi=dpi)
    else
        throw(ArgumentError("context only accepts vector."))
    end
    savefig(p, joinpath(savepath, "pca_$(title)"), save)
    return p
end

function umap(umap1::AbstractVector, umap2::AbstractVector;
              savepath::String=@__DIR__, title="", save=[:svg, :png],
              figsize=(1000, 600), dpi=300,
    )
    p = scatter(umap1, umap2; group=:cell, xlabel="UMAP1", ylabel="UMAP2", color_palette=:glasbey_hv_n256,
        legend=:outertopright, markersize=2, markerstrokewidth=0,
        size=figsize, dpi=dpi, thickness_scaling=2, widen=false)
    savefig(p, joinpath(savepath, "umap_$(title)"), save)
    return p
end
