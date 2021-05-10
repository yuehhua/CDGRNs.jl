SSE(::NullRegression, X::AbstractVecOrMat, y::AbstractVector{T}) where {T} = typemax(float(T))
SSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = residual(model, X, y)'*residual(model, X, y)


MSE(::NullRegression, X::AbstractVecOrMat, y::AbstractVector{T}) where {T} = typemax(float(T))
function MSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector; corrected=false)
    sse = SSE(model, X, y)
    return corrected ? sse / dof(model, kind=:residual) : sse / nobs(model)
end

function MSE(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector{T}) where {K,T}
    clst = predict_cluster(model, X, y)
    total_sse = zero(T)
    for k = 1:K
        m = model.models[k]
        total_sse += SSE(m, _view(X, clst .== k), view(y, clst .== k))
    end
    return total_sse / length(y)
end


likelihood(::NullRegression, X::AbstractVecOrMat, y::AbstractVector{T}) where {T} = zero(float(T))
function likelihood(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector)
    normal = Normal(0, std(model))
    return pdf.(normal, residual(model, X, y))
end


"""
    loglikelihood(model)
    loglikelihood(model, X, y)

Log likelihood of a model.

Returns log likelihood values. If X, y are not given, log likelihood is evaluated by training data.
"""
loglikelihood(::NullRegression; average::Bool=false) = -Inf
loglikelihood(::NullRegression, X::AbstractVecOrMat, y::AbstractVector{T}; average::Bool=false) where {T} = typemin(float(T))
function loglikelihood(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector; average::Bool=false)
    normal = Normal(0, std(model))
    lls = logpdf.(normal, residual(model, X, y))
    return average ? mean(lls) : sum(lls)
end

loglikelihood(model::MixtureRegression{K}; average::Bool=false) where {K} = _loglikelihood(model.likelihoods, model.clusters, K, average)

function loglikelihood(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector; average::Bool=false) where {K}
    likelihoods = _likelihood(model, X, y)
    clusters = map(hard_split, likelihoods...)
    return _loglikelihood(likelihoods, clusters, K, average)
end

function _loglikelihood(ls, clusters, K, avg)
    ll = if K == 1
        sum(log.(ls[1]))
    else
        sum(k -> sum(log.(ls[k][clusters .== k])), 1:K)
    end
    return avg ? ll/length(clusters) : ll
end


"""
    nll(model, X, y)

Negative log likelihood for a model.

Returns negative log likelihood for a linear regression model, which are evaluated by data.
"""
nll(model::AbstractRegression, X::AbstractVecOrMat, y::AbstractVector) = -loglikelihood(model, X, y)


"""
    aic(model)
    aic(model, X, y)

Akaike information criterion of a mixture model.

Returns Akaike information criterion values. If X, y are not given, AIC is evaluated by training data.
"""
aic(model::AbstractRegression; λ=2e-2) = 2(ncoef(model) - λ*loglikelihood(model))
function aic(model::AbstractRegression, X::AbstractVecOrMat, y::AbstractVector; λ=2e-2, kind=:likelihood)
    k = ncoef(model)
    n = length(y)
    if kind == :likelihood
        return 2(k - λ*loglikelihood(model, X, y))
    elseif kind == :mse
        return 2k + n*log(MSE(model, X, y))
    else
        throw(ArgumentError("only :likelihood and :mse available for `kind`."))
    end
end


"""
    bic(model)
    bic(model, X, y)

Bayesian information criterion of a mixture model.

Returns Bayesian information criterion values. If X, y are not given, BIC is evaluated by training data.
"""
bic(model::AbstractRegression; λ=2e-2) = ncoef(model)*log(length(y)) - 2λ*loglikelihood(model)
function bic(model::AbstractRegression, X::AbstractVecOrMat, y::AbstractVector; λ=2e-2, kind=:likelihood)
    k = ncoef(model)
    n = length(y)
    if kind == :likelihood
        return k*log(n) - 2λ*loglikelihood(model, X, y)
    elseif kind == :mse
        return k*log(n) + n*log(MSE(model, X, y))
    else
        throw(ArgumentError("only :likelihood and :mse available for `kind`."))
    end
end
