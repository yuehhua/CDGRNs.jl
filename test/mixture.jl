using DataFrames
using GLM

@testset "mixture" begin
    n = 1000
    k = 2
    iter = 5

    w₁ = 3.
    w₂ = -3.
    b = 1.

    x₁ = randn(n) * 10
    x₂ = randn(n) * 10

    y₁ = w₁.*x₁ .+ b .+ 5*randn(n)
    y₂ = w₂.*x₂ .+ b .+ 5*randn(n)

    data1 = DataFrame(X=x₁, Y=y₁)
    model1_true = lm(@formula(Y ~ X), data1)
    data2 = DataFrame(X=x₂, Y=y₂)
    model2_true = lm(@formula(Y ~ X), data2)

    data = DataFrame(X=vcat(x₁, x₂), Y=vcat(y₁, y₂))
    data.cluster = rand([1,2], 2*n)

    model = fit(MixtureRegression{2}, Matrix(data.X'), data.Y; max_iter=iter, init_clusters=data.cluster)
    @test check_confint(model.models[1].β, stderror(model.models[1]), GLM.coef(model1_true)) | 
        check_confint(model.models[1].β, stderror(model.models[1]), GLM.coef(model2_true))
    @test check_confint(model.models[2].β, stderror(model.models[2]), GLM.coef(model1_true)) | 
        check_confint(model.models[2].β, stderror(model.models[2]), GLM.coef(model2_true))
end
