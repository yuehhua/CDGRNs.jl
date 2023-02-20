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
