function build_tree(prof, true_reg_pairs; linkage=:ward, celltype=nothing, time=nothing)
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
    return tree, D, features
end

function extract_context!(cell_clusters, tree, k)
    cell_clusters[:, Symbol("k$k")] = cutree(tree; k=k)
    return cell_clusters
end
