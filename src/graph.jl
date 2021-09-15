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
