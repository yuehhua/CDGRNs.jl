using Flux: mse

d = 10
numGene = 20
numCell = 30


A = adjacency_list(rand([0,1], numGene, numGene))
train_X = rand(Float32, d, numGene, numCell)

model = Chain(Concentration(d=>numGene),
              GeneRegulatory(A, numGene))

loss(X, y) = mse(model(X), y)

# test model
@show Y = model(train_X)

# test loss
# @show loss(train_X, train_y)

# test gradient
# @show gradient(X -> loss(X, train_y), train_X)