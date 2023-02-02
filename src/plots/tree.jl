function clustermap(D::AbstractMatrix, filename::String; method::String="ward", labels=nothing,
    figsize=(12, 9), cmap="YlGnBu", col_palette=:default, col_colors=nothing, dpi=300)
    fastcluster = pyimport("fastcluster")

    link = fastcluster.linkage(D, method=method)

    if isnothing(labels)
        clustermap(D, link, filename, figsize=figsize, cmap=cmap, dpi=dpi)
    else
        isnothing(col_colors) && (col_colors = generate_plots_colors(labels, col_palette=col_palette))
        clustermap(D, link, col_colors, filename, figsize=figsize, cmap=cmap, dpi=dpi)
    end
end

function clustermap(D::AbstractMatrix, link, col_colors, filename::String;
    figsize=(12, 9), cmap="YlGnBu", dpi=300)
    sns = pyimport("seaborn")
    fig = sns.clustermap(D, row_linkage=link, col_linkage=link, col_colors=col_colors,
            figsize=figsize, cmap=cmap)
    fig.savefig(filename, dpi=dpi)
end

function clustermap(D::AbstractMatrix, link, filename::String;
    figsize=(12, 9), cmap="YlGnBu", dpi=300)
    sns = pyimport("seaborn")
    fig = sns.clustermap(D, row_linkage=link, col_linkage=link, figsize=figsize, cmap=cmap)
    fig.savefig(filename, dpi=dpi)
end

function dendrogram(D::AbstractMatrix, linkage=:ward)
    hc_col = hclust(D, linkage=linkage, branchorder=:optimal)
    return plot(hc_col, xticks=false, yticks=false)
end
