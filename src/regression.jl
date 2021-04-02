using LinearAlgebra
using Distributions
using Statistics

import Statistics: std

mutable struct LinearRegression{V<:AbstractVector,T<:Real}
    β::V
    σ::T
end

LinearRegression(d::Integer) = LinearRegression(rand(d+1), rand())

std(model::LinearRegression) = model.σ

# X ∈ (feat * obs)
design_matrix(x::AbstractVector{T}) where {T} = vcat(ones(T, 1, length(x)), x')
design_matrix(X::AbstractMatrix{T}) where {T} = vcat(ones(T, 1, size(X, 2)), X)

predict(model::LinearRegression, X::AbstractVecOrMat) = design_matrix(X)'*model.β

residual(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = y - predict(model, X)

SSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = residual(model, X, y)'*residual(model, X, y)
MSE(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector) = SSE(model, X, y) / (length(y) - 2)

function fit!(model::LinearRegression, X::AbstractVecOrMat, y::AbstractVector)
    D = design_matrix(X)
    model.β = inv(D*D') * D*y
    r = residual(model, X, y)
    model.σ = sqrt(r'*r / length(r))
    model
end

function fit(::Type{LinearRegression}, x::AbstractVector, y::AbstractVector)
    model = LinearRegression(1)
    fit!(model, x, y)
end

function fit(::Type{LinearRegression}, X::AbstractMatrix, y::AbstractVector)
    model = LinearRegression(size(X, 1))
    fit!(model, X, y)
end

function likelihood(model::LinearRegression, X::AbstractMatrix, y::AbstractVector)
    normal = Normal(0, std(model))
    return pdf.(normal, residual(model, X, y))
end
