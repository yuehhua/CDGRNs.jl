function generate_column_colors(xs::AbstractVector; palette="Set1")
    sns = pyimport("seaborn")
    keys = unique(xs)
    lut = Dict(zip(keys, sns.mpl_palette(palette, length(keys))))
    return map(x -> lut[x], xs)
end

function clustermap(D::AbstractMatrix, labels; method::String="ward", filename="clustermap", ext=".png",
                    figsize=(12, 9), cmap="YlGnBu", dpi=300)
    fastcluster = pyimport("fastcluster")
    sns = pyimport("seaborn")

    link = fastcluster.linkage(D, method=method)
    col_colors = generate_column_colors(labels)
    
    fig = sns.clustermap(D, row_linkage=link, col_linkage=link, col_colors=col_colors,
                         figsize=figsize, cmap=cmap)
    filepath = joinpath(GRN.PROJECT_PATH, "pics", "tf-gene gmm model", "clustering", filename*ext)
    fig.savefig(filepath, dpi=dpi)
end

function dendrogram(D::AbstractMatrix, linkage=:ward)
    hc_col = hclust(D, linkage=linkage, branchorder=:optimal)
    return plot(hc_col, xticks=false, yticks=false)
end
