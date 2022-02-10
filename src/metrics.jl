# SSE(model::NullRegression, X::AbstractVecOrMat, y::AbstractVector) = residual(model, X, y)'*residual(model, X, y)
# SSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = residual(model, X, y)'*residual(model, X, y)


# MSE(::NullRegression, X::AbstractVecOrMat, y::AbstractVector) = var(y)
# function MSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector; corrected=false)
#     sse = SSE(model, X, y)
#     return corrected ? sse / dof(model, kind=:residual) : sse / nobs(model)
# end

# function MSE(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector{T}) where {K,T}
#     clst = predict_cluster(model, X, y)
#     total_sse = zero(T)
#     for k = 1:K
#         m = model.models[k]
#         total_sse += SSE(m, _view(X, clst .== k), view(y, clst .== k))
#     end
#     return total_sse / length(y)
# end


# function likelihood(model::NullRegression, X::AbstractVecOrMat, y::AbstractVector)
#     normal = Normal(0, model.σ)
#     return pdf.(normal, residual(model, X, y))
# end

# function likelihood(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector)
#     normal = Normal(0, std(model))
#     return pdf.(normal, residual(model, X, y))
# end

likelihood(model::AbstractGMR, x::AbstractVector) = pdf.(model.dist, x)
likelihood(model::AbstractGMR, X::AbstractMatrix) = [pdf(model.dist, vec(X[i,:])) for i = 1:size(X, 1)]
likelihood(model::FailedGMR, X::AbstractMatrix{T}) where {T} = [zero(T) for i = 1:size(X, 1)]


"""
    membership(model::GMR[, k_component], X)

Return the probability of membership for an observation.
"""
function membership(model::GMR{K}, X::AbstractMatrix) where {K}
    n = size(X, 1)
    return [pdf(model.dist.components[k], vec(X[i,:])) for i = 1:n, k = 1:K]
end

function membership(model::GMR{K}, k_component, X::AbstractMatrix) where {K}
    @assert k_component ≤ K
    n = size(X, 1)
    return [pdf(model.dist.components[k_component], vec(X[i,:])) for i = 1:n]
end


# """
#     loglikelihood(model[, X, y])

# Log likelihood of a model.

# Returns log likelihood values. If X, y are not given, log likelihood is evaluated by training data.
# """
# loglikelihood(::NullRegression; average::Bool=false) = -Inf
# function loglikelihood(model::NullRegression, X::AbstractVecOrMat, y::AbstractVector; average::Bool=false)
#     n = model.n
#     ll = 0.5 * n * log(2π) + n * log(model.σ)
#     ll -= 0.5 * sum(zscore(y).^2) / n
#     return ll
# end

# function loglikelihood(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector; average::Bool=false)
#     normal = Normal(0, std(model))
#     lls = logpdf.(normal, residual(model, X, y))
#     return average ? mean(lls) : sum(lls)
# end

# loglikelihood(model::MixtureRegression{K}) where {K} = _loglikelihood(model.likelihoods, model.clusters, K)

# function loglikelihood(model::MixtureRegression{K}, X::AbstractVecOrMat, y::AbstractVector) where {K}
#     likelihoods = _likelihood(model, X, y)
#     clusters = map(hard_split, likelihoods...)
#     return _loglikelihood(likelihoods, clusters, K)
# end

loglikelihood(model::AbstractGMR, x::AbstractVector) = logpdf.(model.dist, x)
loglikelihood(model::AbstractGMR, X::AbstractMatrix) = [logpdf(model.dist, vec(X[i,:])) for i = 1:size(X, 1)]
loglikelihood(model::FailedGMR, X::AbstractMatrix) = [-Inf for i = 1:size(X, 1)]

function _loglikelihood(ls, clusters, K)
    ll = if K == 1
        sum(log.(ls[1]))
    else
        n = length(clusters)
        # log(w_k)
        logws = map(k -> log(sum(ls[k])) - log(n), 1:K)
        # log P(y | x, θ)
        logliks = map(k -> sum(log.(ls[k][clusters .== k])), 1:K)
        # log(w_k) + log P(y | x, θ)
        log(sum(exp.(logws + logliks)))
    end
    return ll
end

gen_logpdf(mixtures::Vector{<:AbstractMvNormal}, w::AbstractVector, K::Integer=length(w)) =
    (x...) -> log(sum(exp(map(k -> log(w[k]) + logpdf(mixtures[k], [x...]), 1:K))))


"""
    nll(model, X, y)

Negative log likelihood for a model.

Returns negative log likelihood for a linear regression model, which are evaluated by data.
"""
nll(model::AbstractRegression, X::AbstractVecOrMat, y::AbstractVector) = -loglikelihood(model, X, y)


"""
    aic(model[, X, y])

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

aic(model::AbstractGMR{K}, x::AbstractVector; λ=2e-2) where {K} =
    2(K - λ*sum(loglikelihood(model, x)))
aic(model::AbstractGMR{K}, X::AbstractMatrix; λ=2e-2) where {K} =
    2(K - λ*sum(loglikelihood(model, X)))


"""
    bic(model[, X, y])

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

bic(model::AbstractGMR{K}, x::AbstractVector; λ=2e-2) where {K} =
    K*log(length(x)) - 2λ*sum(loglikelihood(model, x))
bic(model::AbstractGMR{K}, X::AbstractMatrix; λ=2e-2) where {K} =
    K*log(size(X, 1)) - 2λ*sum(loglikelihood(model, X))
