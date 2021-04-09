using LinearAlgebra
using Distributions
using Statistics

import Statistics: std

mutable struct LinearRegression{T<:Real,V<:AbstractVector}
    β::V
    σ::T
    se::V
end

LinearRegression(d::Integer) = LinearRegression(rand(d+1), rand(), rand(d+1))

coef(model::LinearRegression) = model.β
std(model::LinearRegression) = model.σ
stderror(model::LinearRegression) = model.se

# X ∈ (obs × feat)
design_matrix(x::AbstractVector{T}) where {T} = hcat(ones(T, length(x)), x)
design_matrix(X::AbstractMatrix{T}) where {T} = hcat(ones(T, size(X, 1)), X)

predict(model::LinearRegression, X::AbstractVecOrMat) = design_matrix(X)*model.β

residual(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = y - predict(model, X)

SSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = residual(model, X, y)'*residual(model, X, y)
MSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = SSE(model, X, y) / (length(y) - 2)

function fit!(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector{T}) where {T}
    D = design_matrix(X)
    model.β = inv(D'*D) * D'*y
    r = residual(model, X, y)
    model.σ = sqrt(r'*r / length(r))
    # se of slopes
    n = length(y)
    X̄ = mean(X, dims=1)
    ϵₓ = X .- X̄
    SS = sum(x -> x.^2, ϵₓ, dims=1)
    mse = MSE(model, X, y)
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

function likelihood(model::LinearRegression, X::AbstractMatrix, y::AbstractVector)
    normal = Normal(0, std(model))
    return pdf.(normal, residual(model, X, y))
end
