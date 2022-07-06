function pca(pc1::AbstractVector, pc2::AbstractVector, pc3::AbstractVector, labels::AbstractVector;
             savepath::String=@__DIR__, title="", save=[:svg, :png])
    p = Plots.scatter(pc1, pc2, pc3, group=labels, markersize=1, markerstrokewidth=0,
                      xlabel="PC1", ylabel="PC2", zlabel="PC3")
    savefig(p, joinpath(savepath, "pca_$(title)"), save)
    return p
end

function pca(data::AbstractMatrix, labels::AbstractVector, context::AbstractVector)
    α = map(x -> x ? 1. : 0.2, context)
    colors = [RGB(x...) for x in generate_seaborn_colors(labels)]
    p = Plots.scatter(pc[:, 1], pc[:, 2], pc[:, 3], group=labels,
                      markercolor=colors, markeralpha=α, markersize=1, markerstrokewidth=0,
                      xlabel="PC1", ylabel="PC2", zlabel="PC3")
    return p
end

function pca(pc1::AbstractVector, pc2::AbstractVector, labels::AbstractVector,
             context::AbstractVector; savepath::String=@__DIR__, title="", save=[:svg, :png],
             xaxis=1, yaxis=2)
    α = map(x -> x ? 1. : 0.2, context)
    colors = [RGB(x...) for x in generate_seaborn_colors(labels)]
    p = Plots.scatter(pc[:, xaxis], pc[:, yaxis], group=labels, markercolor=colors, markeralpha=α,
                      legend=:outerright, markersize=2, markerstrokewidth=0, dpi=300,
                      size=(1000, 600), thickness_scaling=2, widen=false,
                      xlabel="PC$xaxis", ylabel="PC$yaxis")
    savefig(p, joinpath(savepath, "pca_$(title)"), save)
    return p
end

function umap(umap1::AbstractVector, umap2::AbstractVector;
        savepath::String=@__DIR__, title="", save=[:svg, :png],
    )
    p = scatter(umap1, umap2; group=:cell, xlabel="UMAP1", ylabel="UMAP2", color_palette=:glasbey_hv_n256,
        legend=:outertopright, markersize=2, markerstrokewidth=0,
        dpi=300, size=(1000, 600), thickness_scaling=2, widen=false)
    savefig(p, joinpath(savepath, "umap_$(title)"), save)
    return p
end
