function make_mapping(df::DataFrame, p::Pair{Symbol,Symbol})
    Dict(df[i, p.first] => df[i, p.second] for i = 1:nrow(df))
end

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
