function pca(data::AbstractMatrix; dims=3)
    model = fit(PCA, data; maxoutdim=dims)
    return MultivariateStats.transform(model, data)'
end

umap(data::AbstractMatrix, dims::Int=2; n_neighbors=20, min_dist=0.5) =
    UMAP.umap(data, dims, n_neighbors=n_neighbors, min_dist=min_dist)
