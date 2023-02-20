function to_graph(part_cor::DataFrame; src=:tf, dst=:target,
                  nodes=unique!(vcat(part_cor[:,src], part_cor[:,dst])))
    g = SimpleWeightedDiGraph(length(nodes))
    for row in eachrow(part_cor)
        i = findfirst(nodes .== row[src])
        j = findfirst(nodes .== row[dst])
        add_edge!(g, i, j, row.dist)
    end
    return g
end

function network_entropy(g::AbstractGraph)
    d = vec(sum(Graphs.weights(g), dims=1))
    replace!(d, 0=>eps(0.0))
    N = nv(g)
    return sum(log, d) / (N * log(maximum(d)))
end

# function network_entropy(g::AbstractGraph)
#     d = indegree(g)
#     w = vec(sum(LightGraphs.weights(g), dims=1))
#     S = d ./ w .* log.(w)
#     H = S ./ maximum(S)
#     return mean(H)
# end
