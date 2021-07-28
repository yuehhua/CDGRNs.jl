abstract type AbstractRegression end

struct NullRegression{T<:Real} <: AbstractRegression
    μ::T
    σ::T
    n::Integer
end

mutable struct LinearRegression{T<:Real,V<:AbstractVector} <: AbstractRegression
    β::V
    σ::T
    se::V
    n::Integer
end

LinearRegression(d::Integer) = LinearRegression(rand(d+1), rand(), rand(d+1), 0)

coef(model::NullRegression) = [model.μ, model.σ]
coef(model::LinearRegression) = model.β

std(model::NullRegression) = model.σ
std(model::LinearRegression) = model.σ

stderror(model::LinearRegression) = model.se

ncoef(::NullRegression) = 2
ncoef(model::LinearRegression) = length(model.β)

nobs(model::NullRegression) = model.n
nobs(model::LinearRegression) = model.n

function dof(model::AbstractRegression; kind=:regression)
    if kind == :regression
        return ncoef(model) - 1
    elseif kind == :residual
        return nobs(model) - ncoef(model)
    elseif kind == :total
        return nobs(model) - 1
    else
        throw(ArgumentError("only :regression, :residual and :total available for `kind`."))
    end
end

# X ∈ (obs × feat)
design_matrix(x::AbstractVector{T}) where {T} = hcat(ones(T, length(x)), x)
design_matrix(X::AbstractMatrix{T}) where {T} = hcat(ones(T, size(X, 1)), X)

predict(model::NullRegression, x::AbstractVector) = repeat([model.μ], length(x))
predict(model::NullRegression, X::AbstractMatrix) = repeat([model.μ], size(X, 1))
predict(model::LinearRegression, X::AbstractVecOrMat) = design_matrix(X)*model.β

residual(model::NullRegression, X::AbstractVecOrMat, y::AbstractVector) = y - predict(model, X)
residual(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = y - predict(model, X)

function fit!(model::NullRegression, x::AbstractVector, y::AbstractVector)
    model.μ, model.σ = mean_and_std(y, corrected=true)
    model.n = length(y)
    return model
end

function fit!(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector{T}) where {T}
    model.n = length(y)
    D = design_matrix(X)
    inv_D = inv(D'*D)
    model.β = inv_D * D'*y
    mse = MSE(model, X, y, corrected=true)
    model.σ = sqrt(mse * dof(model, kind=:residual) / nobs(model))
    # se of slopes
    n = length(y)
    X̄ = mean(X, dims=1)
    ϵₓ = X .- X̄
    SS = sum(x -> x.^2, ϵₓ, dims=1)
    se_slope = sqrt.(mse ./ vec(SS))
    # se of intercept
    X̄ = [one(T), (X̄.^2)...]
    SS = [n, SS...]
    se_intercept = sqrt(mse * sum(X̄ ./ SS))
    model.se = [se_intercept, se_slope...]
    model
end

fit(::Type{NullRegression}, x::AbstractVector, y::AbstractVector) = NullRegression(mean_and_std(y, corrected=true)..., length(y))

function fit(::Type{LinearRegression}, x::AbstractVector, y::AbstractVector)
    model = LinearRegression(1)
    fit!(model, x, y)
end

function fit(::Type{LinearRegression}, X::AbstractMatrix, y::AbstractVector)
    model = LinearRegression(size(X, 2))
    fit!(model, X, y)
end


## Gaussian mixture regression

abstract type AbstractGMR{K} <: AbstractRegression end

struct NullGMR{D<:AbstractMvNormal} <: AbstractGMR{1}
    dist::D
    n::Integer
end

struct GMR{K,D<:AbstractMixtureModel} <: AbstractGMR{K}
    dist::D
    n::Integer
end

struct FailedGMR{K} <: AbstractGMR{K}
    n::Integer
end

function fit(::Type{GMR{1}}, X::AbstractMatrix)
    μ, σ = mean_and_std(X, 1, corrected=true)
    mvn = MvNormal(vec(μ), diagm(vec(σ)))
    return NullGMR(mvn, size(X, 1))
end

function fit(::Type{GMR{K}}, X::AbstractMatrix) where {K}
    try
        gmm = GaussianMixtures.GMM(K, X, kind=:full)
        model = MixtureModel(gmm)
        return GMR{K,typeof(model)}(model, size(X, 1))
    catch e
        if e isa PosDefException
            return FailedGMR{K}(size(X, 1))
        else
            rethrow(e)
        end
    end
end

assign_clusters(::NullGMR, X::AbstractMatrix) = ones(Int, size(X, 1))

function assign_clusters(model::GMR, X::AbstractMatrix)
    posterior = membership(model, X)
    return assign_clusters(posterior)
end


# population correlation

correlation(σ_xy, σ_xx, σ_yy) = σ_xy / sqrt(σ_xx * σ_yy)

correlation(::FailedGMR) = [0.]

function correlation(model::NullGMR)
    Σ = model.dist.Σ
    [correlation(Σ[1,2], Σ[1,1], Σ[2,2])]
end

correlation(model::GMR) = [correlation(c.Σ[1,2], c.Σ[1,1], c.Σ[2,2]) for c in model.dist.components]


# sample correlation

function correlation(xs::AbstractVector, ys::AbstractVector, cluster::AbstractVector)
    ks = unique(Array(cluster))
    ρs = similar(xs, maximum(ks))
    for clst in ks
        sel = clst .== cluster
        ρs[clst] = cor(xs[sel], ys[sel])
    end
    return ρs
end

fisher_transform(ρs::AbstractVector) = sqrt(length(ρs)-3) .* atanh.(ρs)
