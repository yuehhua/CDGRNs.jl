function generate_column_colors(xs::AbstractVector; palette="rainbow")
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
    filepath = joinpath(CDGRN.PROJECT_PATH, "pics", "tf-gene gmm model", "clustering", filename*ext)
    fig.savefig(filepath, dpi=dpi)
end

function dendrogram(D::AbstractMatrix, linkage=:ward)
    hc_col = hclust(D, linkage=linkage, branchorder=:optimal)
    return plot(hc_col, xticks=false, yticks=false)
end

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
    colors = generate_column_colors(labels)
    colors = map(x -> RGB(x...), colors)
    p = Plots.scatter(pc[:, 1], pc[:, 2], pc[:, 3], group=labels,
                      markercolor=colors, markeralpha=α, markersize=1, markerstrokewidth=0,
                      xlabel="PC1", ylabel="PC2", zlabel="PC3")
    return p
end

function plot_regulations(data::DataFrame, x::String, y::String, z::String; model=nothing, spliced=false)
    plotly()

    x_sym = Symbol(x * "_s")
    y_sym = Symbol(y * "_s")
    z_sym = spliced ? Symbol(z * "_s") : Symbol(z * "_u")
    xlabel = "Spliced mRNA of $x"
    ylabel = "Spliced mRNA of $y"
    zlabel = spliced ? "Spliced mRNA of $z" : "Unspliced mRNA of $z"
    p = scatter(data[:,x_sym], data[:,y_sym], data[:,z_sym], group=data[:,:cell_type],
                markersize=2, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
    if !isnothing(model)
        x = range(floor(minimum(data[:,x_sym])), stop=ceil(maximum(data[:,x_sym])), length=100)
        y = range(floor(minimum(data[:,y_sym])), stop=ceil(maximum(data[:,y_sym])), length=100)
        plot!(p, x, y, (x,y) -> predict(model, [x, y]'), st=:surface,
              c=:grey, fillalpha=0.7, cb=false)
    end
    return p
end

function plot_regulations(data::DataFrame, x::String, y::String; model=nothing, spliced=false)
    plotly()

    x_sym = Symbol(x * "_s")
    y_sym = spliced ? Symbol(y * "_s") : Symbol(y * "_u")
    xlabel = "Spliced mRNA of $x"
    ylabel = spliced ? "Spliced mRNA of $y" : "Unspliced mRNA of $y"
    p = scatter(data[:,x_sym], data[:,y_sym], group=data[:,:cell_type],
                markersize=3, xlabel=xlabel, ylabel=ylabel)
    if !isnothing(model)
        x = range(floor(minimum(data[:,x_sym])), stop=ceil(maximum(data[:,x_sym])), length=100)
        plot!(p, x, y -> predict(model, [y]'), st=:surface,
              c=:grey, fillalpha=0.7, cb=false)
    end
    return p
end

function make_vis_data(target_gene, tfs, prof, pairs; cluster=:)
    tf_set = unique(pairs[pairs.target .== target_gene, :tf])
    data = DataFrame()
    data.cell_type = prof.obs.clusters
    data[:, Symbol(target_gene * "_u")] = vec(get_gene_expr(prof, target_gene, :Mu))
    data[:, Symbol(target_gene * "_s")] = vec(get_gene_expr(prof, target_gene, :Ms))
    for tf in tf_set
        data[:, Symbol(tf * "_s")] = vec(get_gene_expr(tfs, tf, :Ms))
    end
    return data[cluster, :]
end
