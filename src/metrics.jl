SSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = residual(model, X, y)'*residual(model, X, y)


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
function loglikelihood(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector)
    normal = Normal(0, std(model))
    return logpdf.(normal, residual(model, X, y))
end

function loglikelihood(model::MixtureRegression{K}; average::Bool=false) where {K}
    ll = sum(k -> sum(log.(model.likelihoods[k][model.clusters .== k])), 1:K)
    return average ? ll/length(model.clusters) : ll
end

function loglikelihood(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector; average::Bool=false) where {K}
    likelihoods = _likelihood(model, X, y)
    clusters = map(hard_split, likelihoods...)
    ll = sum(k -> sum(log.(likelihoods[k][clusters .== k])), 1:K)
    return average ? ll/length(y) : ll
end


"""
    nll(model, X, y)

Negative log likelihood for a model.

Returns negative log likelihood for a linear regression model, which are evaluated by data.
"""
nll(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = -loglikelihood(model, X, y)


"""
    aic(model)
    aic(model, X, y)

Akaike information criterion of a mixture model.

Returns Akaike information criterion values. If X, y are not given, AIC is evaluated by training data.
"""
aic(model::MixtureRegression; λ=2e-2) = 2(ncoef(model) - λ*loglikelihood(model))
function aic(model::MixtureRegression, X::AbstractVecOrMat, y::AbstractVector; λ=2e-2, kind=:likelihood)
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
bic(model::MixtureRegression; λ=2e-2) = ncoef(model)*log(length(y)) - 2λ*loglikelihood(model)
function bic(model::MixtureRegression, X::AbstractVecOrMat, y::AbstractVector; λ=2e-2, kind=:likelihood)
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
