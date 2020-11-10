using LinearAlgebra
using Flux
using Flux: glorot_uniform, mse
using Flux: @functor
using GeometricFlux
using GraphSignals

d = 10
numGene = 20
numCell = 30

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
    r = [softplus.(view(l.W,:,i)' * view(X,:,i,:) .+ view(l.b,i)) for i = 1:size(l.W,2)]
    vcat(r...)
end



## Gene Regulatory Layer

struct GeneRegulatory{V<:AbstractFeaturedGraph,T<:AbstractVector,S<:AbstractMatrix} <: MessagePassing
    fg::V
    W::T
    b::T
    R::S
end

GeneRegulatory(al::AbstractVector{<:AbstractVector{<:Integer}}, dim::Integer;
               init=glorot_uniform, T::DataType=Float32) =
    GeneRegulatory(FeaturedGraph(al), dim, init=init, T=T)

GeneRegulatory(dim::Integer; init=glorot_uniform, T::DataType=Float32) = GeneRegulatory(NullGraph(), dim, init=init, T=T)

function GeneRegulatory(fg::AbstractFeaturedGraph, dim::Integer;
                        init=glorot_uniform, T::DataType=Float32)
    W = T.(init(dim))
    b = T.(init(dim))
    R = T.(init(dim, dim))
    R .-= LowerTriangular(R)
    R[1,1] = init(1)[1]
    GeneRegulatory(fg, W, b, R)
end

@functor GeneRegulatory

function update_batch_edge(l::GeneRegulatory, adj, X::AbstractMatrix)
    n = size(adj, 1)
    vcat([apply_batch_message(l, i, adj[i], X) for i in 1:n]...)
end

@inline function apply_batch_message(l::GeneRegulatory, i, js, X::AbstractMatrix)
    h = σ.(view(l.W,js) .* log.(view(X,js,:)) .+ view(l.b,js))
    m = [view(h,:,j)'*view(l.R,js,js)*view(h,:,j) for j = 1:size(h,2)]
    m'
end

propagate(l::GeneRegulatory, adj, X::AbstractMatrix) = update_batch_edge(l, adj, X)

function (l::GeneRegulatory)(fg::FeaturedGraph)
    X = propagate(l, adjacency_list(fg), node_feature(fg))
    FeaturedGraph(graph(fg), X)
end

function (l::GeneRegulatory)(X::AbstractMatrix)
    @assert has_graph(l.fg) "A GeneRegulatory created without a graph must be given a FeaturedGraph as an input."
    propagate(l, adjacency_list(l.fg), X)
end

(l::GeneRegulatory)(adj, X::AbstractMatrix) = propagate(l, adj, X)


A = adjacency_list(rand([0,1], numGene, numGene))
train_X = rand(Float32, d, numGene, numCell)

model = Chain(Concentration(d=>numGene),
              GeneRegulatory(A, numGene))

loss(X, y) = mse(model(X), y)

# test model
@show Y = model(train_X)

# test loss
# @show loss(train_X, train_y)

# test gradient
# @show gradient(X -> loss(X, train_y), train_X)