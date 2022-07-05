function build_tree(prof, true_reg_pairs; linkage=:ward, celltype=:clusters,
                    time=:latent_time, col_palette=:default, save=nothing)
    features = DataFrame()
    isnothing(celltype) || (features[!, :cell] = prof.obs[!, celltype])
    isnothing(time) || (features[!, :time] = prof.obs[!, time])
    for res in true_reg_pairs
        colname = res[:tf_name] * "_" * res[:gene_name]
        features[!, colname] = res[:clusters]
    end
    
    data = Array(features[:, 3:end])
    D = pairwise(Hamming(), data, dims=1)
    tree = hclust(D, linkage=linkage, branchorder=:optimal)
    isnothing(save) || clustermap(D, features.cell, save, col_palette=col_palette)
    return tree, features
end

function extract_context!(cell_clusters, tree, k)
    cell_clusters[:, Symbol("k$k")] = cutree(tree; k=k)
    return cell_clusters
end
