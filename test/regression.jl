using GLM
using Statistics

T = Float64
n = 100
d = 5

@testset "regression" begin
    # simple regression
    w = 3.0
    b = -1.0
    σ = 2.0

    x = T.(collect(1:n))
    y = w.*x .+ b .+ randn(n)*σ
    df = DataFrame(X=x, Y=y)

    @test size(design_matrix(x)) == (n, 2)
    @test design_matrix(x)[:, 1] == ones(T, length(x))

    model = GRN.fit(LinearRegression, x, y)
    model_true = lm(@formula(Y ~ X), df)
    @test size(model.β) == (2,)
    @test model.β ≈ GLM.coef(model_true)
    @test check_confint(model.β, stderror(model), GLM.coef(model_true))
    @test std(model) ≈ Statistics.std(residual(model, x, y), corrected=false)

    @test size(predict(model, x)) == size(y)
    @test size(residual(model, x, y)) == size(y)


    # multiple regression
    w = [1., 2., 3., 4., 5]

    X = randn(n, d)
    y = X*w .+ b .+ randn(n)*σ
    df = DataFrame(X1=X[:,1], X2=X[:,2], X3=X[:,3], X4=X[:,4], X5=X[:,5], Y=y)

    @test size(design_matrix(X)) == (n, d+1)
    @test design_matrix(X)[:, 1] == ones(T, size(X, 1))

    model = fit(LinearRegression, X, y)
    model_true = lm(@formula(Y ~ X1 + X2 + X3 + X4 + X5), df)
    @test size(model.β) == (d+1,)
    @test model.β ≈ GLM.coef(model_true)
    @test check_confint(model.β, stderror(model), GLM.coef(model_true))
    @test std(model) ≈ Statistics.std(residual(model, X, y), corrected=false)

    @test size(predict(model, X)) == size(y)
    @test size(residual(model, X, y)) == size(y)
end

@testset "GMR" begin
    # Gaussian mixture regression

    w = [1., 2., 3., 4., 5]
    b = -1.0
    σ = 2.0

    X = randn(n, d)
    y = X*w .+ b .+ randn(n)*σ
    X = hcat(X, y)

    nullgmr = fit(GMR{1}, X)
    @test nullgmr isa NullGMR
    @test typeof(nullgmr) <: AbstractGMR{1}

    gmr = fit(GMR{3}, X)
    @test gmr isa GMR
    @test typeof(gmr) <: AbstractGMR{3}
end
