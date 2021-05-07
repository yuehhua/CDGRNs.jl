using GLM
using GaussianMixtures

struct MixtureRegression{K,R,T<:Integer,S}
    models::Vector{R}
    clusters::Vector{T}
    likelihoods::Vector{S}

    function MixtureRegression{K}(models::Vector{R}, clusters::Vector{T}, likelihoods::AbstractVector) where {K,R,T}
        new{K,R,T,eltype(likelihoods)}(models, clusters, likelihoods)
    end
end

ncoef(model::MixtureRegression) = sum(ncoef, model.models)

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

_view(X::AbstractMatrix, idx) = view(X, idx, :)
_view(x::AbstractVector, idx) = view(x, idx)

function maximize_likelihood!(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector) where {K}
    for c = 1:K
        X_c = _view(X, model.clusters .== c)
        y_c = view(y, model.clusters .== c)
        fit!(model.models[c], X_c, y_c)
        model.likelihoods[c] .= likelihood(model.models[c], X, y)
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

function fit!(model::MixtureRegression, X::AbstractVecOrMat, y::AbstractVector{T};
              max_iter::Integer=5) where {T<:Real}
    for i in 1:max_iter
        maximize_likelihood!(model, X, y)
        update_expectation!(model)
    end
    maximize_likelihood!(model, X, y)
    update_expectation!(model; methods=hard_split)
    return model
end

function fit(::Type{MixtureRegression{K}}, X::AbstractVecOrMat, y::AbstractVector{T};
             max_iter::Integer=5, init=()->gmm_init(K, X, y)) where {K,T<:Real}
    # if K == 1
    #     models = [fit(LinearRegression, X, y)]
    #     clusters = ones(T, size(y))
    #     likelihoods = [likelihood(models[1], X, y)]
    #     model = MixtureRegression{K}(models, clusters, likelihoods)
    #     return model
    # else
        n = length(y)
        models = [LinearRegression(size(X, 2)) for i = 1:K]
        init_clusters = init()
        init_likelihoods = [Vector{T}(undef, n) for i = 1:K]
        model = MixtureRegression{K}(models, init_clusters, init_likelihoods)
        return fit!(model, X, y; max_iter=max_iter)
    # end
end

random_init(k, n) = rand(collect(1:k), n)

function gmm_init(k::Integer, xs::AbstractVector, y::AbstractVector)
    train = hcat(xs, y)
	gmm = GaussianMixtures.GMM(k, train, kind=:full)
	ll = GaussianMixtures.llpg(gmm, train)
	clusters = vec(map(x -> x[2], argmax(ll, dims=2)))
    return clusters
end

function gmm_init(k::Integer, X::AbstractMatrix, y::AbstractVector)
    train = hcat(X', y)
	gmm = GaussianMixtures.GMM(k, train, kind=:full)
	ll = GaussianMixtures.llpg(gmm, train)
	clusters = vec(map(x -> x[2], argmax(ll, dims=2)))
    return clusters
end

predict(model::MixtureRegression{K}, X::AbstractMatrix) where {K} = [predict(model.models[k], X) for k = 1:K]

_likelihood(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector) where {K} = [likelihood(model.models[k], X, y) for k = 1:K]

predict_cluster(model::MixtureRegression, X::AbstractVecOrMat, y::AbstractVector) = map(hard_split, _likelihood(model, X, y)...)
