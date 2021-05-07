using Statistics

import Statistics: std

mutable struct LinearRegression{T<:Real,V<:AbstractVector}
    β::V
    σ::T
    se::V
    n::Integer
end

LinearRegression(d::Integer) = LinearRegression(rand(d+1), rand(), rand(d+1), 0)

coef(model::LinearRegression) = model.β
std(model::LinearRegression) = model.σ
stderror(model::LinearRegression) = model.se
ncoef(model::LinearRegression) = length(model.β)
nobs(model::LinearRegression) = model.n

function dof(model::LinearRegression; kind=:regression)
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

predict(model::LinearRegression, X::AbstractVecOrMat) = design_matrix(X)*model.β

residual(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = y - predict(model, X)

function fit!(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector{T}) where {T}
    model.n = length(y)
    D = design_matrix(X)
    model.β = inv(D'*D) * D'*y
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

function fit(::Type{LinearRegression}, x::AbstractVector, y::AbstractVector)
    model = LinearRegression(1)
    fit!(model, x, y)
end

function fit(::Type{LinearRegression}, X::AbstractMatrix, y::AbstractVector)
    model = LinearRegression(size(X, 2))
    fit!(model, X, y)
end
