## Protein Concentration Layer

struct Concentration{T<:AbstractMatrix,S<:AbstractVector}
    W::T
    b::S
end

function Concentration(ch::Pair{<:Integer,<:Integer}; init=glorot_uniform, T::DataType=Float32)
    W = T.(init(ch[1], ch[2]))
    b = T.(init(ch[2]))
    Concentration(W, b)
end

@functor Concentration

function (l::Concentration)(X::AbstractArray)
    WX = reshape(sum(l.W .* X, dims=1), size(X,2), size(X,3))
    softplus.(WX .+ l.b)
end



## Gene Regulatory Layer

struct GeneRegulatory{V<:AbstractFeaturedGraph,T<:AbstractVector} <: MessagePassing
    fg::V
    W::T
    b::T
    α::T
    β::T
end

GeneRegulatory(al::AbstractVector{<:AbstractVector{<:Integer}}, dim::Integer;
               init=glorot_uniform, T::DataType=Float32) =
    GeneRegulatory(FeaturedGraph(al), dim, init=init, T=T)

GeneRegulatory(dim::Integer; init=glorot_uniform, T::DataType=Float32) = GeneRegulatory(NullGraph(), dim, init=init, T=T)

function GeneRegulatory(fg::AbstractFeaturedGraph, dim::Integer;
                        init=glorot_uniform, T::DataType=Float32)
    W = T.(init(dim))
    b = T.(init(dim))
    α = T.(init(dim))
    β = T.(init(dim))
    GeneRegulatory(fg, W, b, α, β)
end

@functor GeneRegulatory

function update_batch_edge(l::GeneRegulatory, adj, X::AbstractMatrix)
    H = @. σ(l.W * log(X) + l.b)
    n = size(adj, 1)
    vcat([apply_batch_message(l, adj[i], H) for i in 1:n]...)
end

@inline function apply_batch_message(l::GeneRegulatory, js, H::AbstractMatrix)
    M = view(l.β,js) .+ view(l.α,js) .* view(H,js,:)
    prod(M, dims=1)
end

propagate(l::GeneRegulatory, adj, X::AbstractMatrix) = update_batch_edge(l, adj, X)

function (l::GeneRegulatory)(fg::FeaturedGraph)
    X = propagate(l, adjacency_list(fg), node_feature(fg))
    FeaturedGraph(graph(fg), nf=X)
end

function (l::GeneRegulatory)(X::AbstractMatrix)
    @assert has_graph(l.fg) "A GeneRegulatory created without a graph must be given a FeaturedGraph as an input."
    propagate(l, adjacency_list(l.fg), X)
end

(l::GeneRegulatory)(adj, X::AbstractMatrix) = propagate(l, adj, X)