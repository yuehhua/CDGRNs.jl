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

make_pairset(reg::DataFrame) = Set(convert_pairs(uppercase, reg.tf, reg.target))

function query_pairset(corr_pairs::DataFrame, reg_pairs::Set)
    query_pairs = convert_pairs(uppercase, corr_pairs.tf, corr_pairs.target)
    return [x in reg_pairs for x in query_pairs]
end

convert_pairs(f, xs, ys) = [(f(x), f(y)) for (x, y) in zip(xs, ys)]

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
