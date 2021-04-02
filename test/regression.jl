using GLM
using Statistics

@testset "regression" begin
    T = Float64
    n = 100
    d = 5

    # simple regression
    w = 3.0
    b = -1.0
    σ = 2.0

    x = T.(collect(1:n))
    y = w.*x .+ b .+ randn(n)*σ
    df = DataFrame(X=x, Y=y)

    @test size(design_matrix(x)) == (2, n)
    @test design_matrix(x)[1, :] == ones(T, length(x))

    model = GRN.fit(LinearRegression, x, y)
    model_true = lm(@formula(Y ~ X), df)
    @test size(model.β) == (2,)
    @test model.β ≈ GLM.coef(model_true)
    @test std(model) ≈ Statistics.std(residual(model, x, y), corrected=false)

    @test size(predict(model, x)) == size(y)
    @test size(residual(model, x, y)) == size(y)


    # multiple regression
    w = [1., 2., 3., 4., 5]

    X = randn(d, n)
    y = X'*w .+ b .+ randn(n)*σ
    df = DataFrame(X1=X[1,:], X2=X[2,:], X3=X[3,:], X4=X[4,:], X5=X[5,:], Y=y)

    @test size(design_matrix(X)) == (d+1, n)
    @test design_matrix(X)[1, :] == ones(T, size(X, 2))

    model = fit(LinearRegression, X, y)
    model_true = lm(@formula(Y ~ X1 + X2 + X3 + X4 + X5), df)
    @test size(model.β) == (d+1,)
    @test model.β ≈ GLM.coef(model_true)
    @test std(model) ≈ Statistics.std(residual(model, X, y), corrected=false)

    @test size(predict(model, X)) == size(y)
    @test size(residual(model, X, y)) == size(y)
end