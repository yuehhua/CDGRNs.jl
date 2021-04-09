using Distributions
using GLM
using StatsBase

struct MixtureRegression{K,R,T<:Integer,S}
    models::Vector{R}
    clusters::Vector{T}
    likelihoods::Vector{S}
end

hard_split(likelihoods...) = argmax.(zip(likelihoods...))

function probabilistic_split(likelihoods...)
    xs = collect(1:length(likelihoods))
    res = []
    for w in zip(likelihoods...)
        w = collect(w ./ sum(w))
        selected = any(w .> 0.9) ? xs[w .> 0.9][1] : StatsBase.sample(xs, Weights(collect(w)))
        push!(res, selected)
    end
    res
end

function maximize_likelihood!(model::MixtureRegression{K}, X::AbstractMatrix, y::AbstractVector) where {K}
    if isempty(model.models)
        for c = 1:K
            X_c = X[:, model.clusters .== c]
            y_c = y[model.clusters .== c]
            push!(model.models, fit(LinearRegression, X_c, y_c))
            push!(model.likelihoods, likelihood(model.models[c], X, y))
        end
    else
        for c = 1:K
            X_c = X[:, model.clusters .== c]
            y_c = y[model.clusters .== c]
            fit!(model.models[c], X_c, y_c)
            model.likelihoods[c] .= likelihood(model.models[c], X, y)
        end
    end
    return model.likelihoods
end

function update_expectation!(model::MixtureRegression; methods=probabilistic_split)
    if typeof(methods) == typeof(probabilistic_split)
        model.clusters .= probabilistic_split(model.likelihoods...)
        return map(x -> x[1], model.clusters)
    elseif typeof(methods) == typeof(hard_split)
        return map(hard_split, model.likelihoods...)
    end
end

function fit!(model::MixtureRegression, X::AbstractMatrix, y::AbstractVector{T};
              max_iter::Integer=5) where {T<:Real}
    for i in 1:max_iter
        maximize_likelihood!(model, X, y)
        update_expectation!(model)
    end
    maximize_likelihood!(model, X, y)
    update_expectation!(model; methods=hard_split)
    return model
end

function fit(::Type{MixtureRegression{K}}, X::AbstractMatrix, y::AbstractVector{T};
             max_iter::Integer=5, init=random_init) where {K,T<:Real}
    S = Vector{T}
    init_clusters = init(K, length(y))
    model = MixtureRegression{K,LinearRegression,Int,S}(LinearRegression[], init_clusters, S[])
    return fit!(model, X, y; max_iter=max_iter)
end

random_init(k, n) = rand(collect(1:k), n)
