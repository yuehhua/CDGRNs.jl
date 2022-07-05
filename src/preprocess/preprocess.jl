function get_regulation_expr(prof::Profile, tfs::Profile, true_reg::DataFrame; labels=:clusters)
    tf_list = unique(true_reg.tf)
    gene_list = unique(true_reg.target)
    
    df = DataFrame(cell=prof.obs[!,labels], time=prof.obs.latent_time)
    for gene in gene_list
        df[:, Symbol(gene * "_u")] = vec(get_gene_expr(prof, gene, :Mu))
    end
    for tf in tf_list
        df[:, Symbol(tf * "_s")] = vec(get_gene_expr(tfs, tf, :Ms))
    end
    return df
end

function build_graph(df::DataFrame)
    gene_set = unique!(vcat(df.tf, df.target))
    N = length(gene_set)
    sg = SimpleDiGraph(N)
    for (tf, targ) in zip(df.tf, df.target)
        i = findfirst(gene_set .== tf)
        j = findfirst(gene_set .== targ)
        add_edge!(sg, i, j)
    end
    return sg
end
