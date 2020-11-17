function make_graph(data::DataFrame, encoding::Dict, graph_size::Int, n::Int=nrow(data))
    dg = SimpleDiGraph{UInt32}(graph_size)
    for i = 1:n
        g = encoding[data[i, :gene_id]]
        tf = encoding[data[i, :tf_id]]
        add_edge!(dg, tf, g)
    end
    dg
end