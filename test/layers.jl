T = Float32
D = 10
numGene = 1000
numCell = 300

A = rand([0,1], numGene, numGene)
A = T.(A .| A')

@testset "layers" begin
    @testset "Concentration" begin
        X = rand(T, D, numGene, numCell)
        l = Concentration(D=>numGene)
        @test size(l.W) == (D, numGene)
        @test size(l.b) == (numGene,)

        Y = l(X)
        @test size(Y) == (numGene, numCell)

        g = Zygote.gradient(x -> sum(l(x)), X)[1]
        @test size(g) == size(X)

        g = Zygote.gradient(model -> sum(model(X)), l)[1]
        @test size(g.W) == size(l.W)
        @test size(g.b) == size(l.b)
    end

    @testset "GeneRegulatory" begin
        X = rand(T, numGene, numCell)
        l = GeneRegulatory(adjacency_list(A), numGene)
        @test size(l.W) == (numGene,)
        @test size(l.b) == (numGene,)
        @test size(l.α) == (numGene,)
        @test size(l.β) == (numGene,)

        Y = l(X)
        @test size(Y) == (numGene, numCell)

        g = Zygote.gradient(x -> sum(l(x)), X)[1]
        @test size(g) == size(X)

        g = Zygote.gradient(model -> sum(model(X)), l)[1]
        @test size(g.W) == size(l.W)
        @test size(g.b) == size(l.b)
        @test size(g.α) == size(l.α)
        @test size(g.β) == size(l.β)
    end
end
