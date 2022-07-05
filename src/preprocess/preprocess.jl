function corr_table(pairs)
    # unit: component
    cor_pairs = map(pairs) do x
        (tf=x[:tf_name],
         target=x[:gene_name],
         best_k=x[:best_k],
         corr=GaussianMixtureRegressions.correlation(x[:model]))
    end
    df = DataFrames.flatten(DataFrame(cor_pairs), :corr)
    df.adjusted_corr = GaussianMixtureRegressions.fisher_transform(df.corr)
    return df
end

function make_pairset(reg::DataFrame)
    return Set(map((x,y) -> (uppercase(x), uppercase(y)), reg.tf, reg.target))
end

function query_pairset(corr_pairs::DataFrame, reg_pairs::Set)
    query_pairs = map((x,y) -> (uppercase(x), uppercase(y)), corr_pairs.tf, corr_pairs.target)
    return map(x -> x in reg_pairs, query_pairs)
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
