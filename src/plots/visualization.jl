function pca_transform(data::AbstractMatrix; dims=3)
    model = fit(PCA, data; maxoutdim=dims)
    return MultivariateStats.transform(model, data)'
end

function plot_3d_pca(data::AbstractMatrix, labels::AbstractVector)
    pc = pca_transform(data)
    p = Plots.scatter(pc[:, 1], pc[:, 2], pc[:, 3], group=labels, markersize=1, markerstrokewidth=0,
                      xlabel="PC1", ylabel="PC2", zlabel="PC3")
    return p
end

function plot_3d_pca(data::AbstractMatrix, labels::AbstractVector, context::AbstractVector)
    pc = pca_transform(data)
    α = map(x -> x ? 1. : 0.2, context)
    colors = generate_seaborn_colors(labels)
    colors = map(x -> RGB(x...), colors)
    p = Plots.scatter(pc[:, 1], pc[:, 2], pc[:, 3], group=labels,
                      markercolor=colors, markeralpha=α, markersize=1, markerstrokewidth=0,
                      xlabel="PC1", ylabel="PC2", zlabel="PC3")
    return p
end

function plot_2d_pca(data::AbstractMatrix, labels::AbstractVector, context::AbstractVector; xaxis=1, yaxis=2)
    pc = pca_transform(data)
    α = map(x -> x ? 1. : 0.2, context)
    colors = generate_seaborn_colors(labels)
    colors = map(x -> RGB(x...), colors)
    p = Plots.scatter(pc[:, xaxis], pc[:, yaxis], group=labels, markercolor=colors, markeralpha=α,
                      legend=:outerright, markersize=2, markerstrokewidth=0, dpi=300,
                      xlabel="PC$xaxis", ylabel="PC$yaxis")
    return p
end
