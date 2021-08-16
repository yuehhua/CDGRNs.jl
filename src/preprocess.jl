function make_graph(data::DataFrame, encoding::Dict, graph_size::Int, n::Int=nrow(data);
                    rev::Bool=false)
    dg = SimpleDiGraph{UInt32}(graph_size)
    for i = 1:n
        g = encoding[data[i, :gene_symbol]]
        tf = encoding[data[i, :tf_symbol]]
        rev ? add_edge!(dg, tf, g) : add_edge!(dg, g, tf)
    end
    dg
end

function get_regulation_expr(prof::Profile, tfs::Profile, true_reg::DataFrame)
    tf_list = unique(true_reg.tf)
    gene_list = unique(true_reg.target)
    
    df = DataFrame(cell=prof.obs.clusters, time=prof.obs.latent_time)
    for gene in gene_list
        df[:, Symbol(gene * "_u")] = vec(get_gene_expr(prof, gene, :Mu))
    end
    for tf in tf_list
        df[:, Symbol(tf * "_s")] = vec(get_gene_expr(tfs, tf, :Ms))
    end
    return df
end
