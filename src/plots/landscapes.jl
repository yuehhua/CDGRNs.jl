logpdf(model::GMR) = (x,y) -> logpdf(model.dist, [x,y])

function likelihood_landscape(df::DataFrame, tf_name, gene_name, model::GMR;
        xlabel="log spliced RNA of TF gene, $(tf_name)",
        ylabel="log unspliced RNA of target gene, $(gene_name)",
        savepath::String=@__DIR__, title="likelihood_landscape", save=[:svg],
        figsize=(800, 600), dpi=300)
    
    xmin, xmax = extrema(df.logX)
    ymin, ymax = extrema(df.logY)
    landscape = logpdf(model)
    p = contour!(xmin:0.01:xmax, ymin:0.01:ymax, landscape, seriescolor=:OrRd)
    scatter!(df.logX, df.logY, xlabel=xlabel, ylabel=ylabel,
        markersize=2, markerstrokewidth=0, legends=false,
        thickness_scaling=2, widen=false, size=figsize, dpi=dpi
    )
    savefig(p, joinpath(savepath, "$(title)_$(tf_name)-$(gene_name)"), save)
    return p
end

function cluster_landscape(df::DataFrame, tf_name, gene_name, clst;
        xlabel="log spliced RNA of TF gene, $(tf_name)",
        ylabel="log unspliced RNA of target gene, $(gene_name)",
        savepath::String=@__DIR__, title="cluster_landscape", save=[:svg],
        figsize=(800, 600), dpi=300)
    
    p = scatter(df.logX, df.logY, color=clst, xlabel=xlabel, ylabel=ylabel,
        markersize=2, markerstrokewidth=0,
        thickness_scaling=2, widen=false, size=figsize, dpi=dpi
    )
    savefig(p, joinpath(savepath, "$(title)_$(tf_name)-$(gene_name)"), save)
    return p
end
