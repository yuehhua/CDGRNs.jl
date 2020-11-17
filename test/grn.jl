using Flux: mse

## Load data

dir = joinpath(@__DIR__, "..", "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)

## Dimensions

d = 10
numGene = 20
numCell = 30

## Preprocessing

A = adjacency_list(rand([0,1], numGene, numGene))
train_X = rand(Float32, d, numGene, numCell)

## Model

model = Chain(Concentration(d=>numGene),
              GeneRegulatory(A, numGene))

## Loss

loss(X, y) = mse(model(X), y)

# test model
@show Y = model(train_X)

# test loss
# @show loss(train_X, train_y)

# test gradient
# @show gradient(X -> loss(X, train_y), train_X)