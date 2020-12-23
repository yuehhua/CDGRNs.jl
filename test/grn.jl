using Flux
using Flux: mse
using GRN
using GraphSignals
using DataFrames
using CSV
using JLD2
using SnowyOwl
using Zygote

## Load data

dir = joinpath(@__DIR__, "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)

gene_list = DataFrame(symbol=uppercase.(prof.var.index))
CSV.write(joinpath(dir, "gene_list.csv"), gene_list)

@load "results/tf_gene_network.jld2" dg
@load "results/gene_set.jld2" gene_set
@load "results/tf_set.jld2" tf_set
@load "results/gene2num.jld2" gene2num

## Dimensions

d = 3
numGene = nv(dg)
numCell = 200#ncol(prof)

## Preprocessing

function truncat_adjl!(adjl::Vector, n)
    adjl = adjl[1:n]
    for i = 1:n
        l = adjl[i]
        adjl[i] = l[l .<= n]
    end
    adjl
end

function truncat_gene2num!(gene2num::Dict, n)
    for (k, v) in gene2num
        if v > n
            delete!(gene2num, k)
        end
    end
    gene2num
end

function prepare_X(prof::Profile, gene2num::Dict, d, numGene, numCell)
    spliced = prof.layers[:spliced]
    velocity = prof.layers[:velocity]
    fit_alpha = prof.var.fit_alpha
    train_X = Array{Float32,3}(undef, d, numGene, numCell)
    for (i, g) = enumerate(prof.var.index), j = 1:numCell
        k = get(gene2num, uppercase(g), 0)
        if k != 0
            α = ismissing(fit_alpha[i]) ? zero(Float32) : fit_alpha[i]
            train_X[:, k, j] .= [spliced[i,j], velocity[i,j], α]
        end
    end
    train_X
end

function prepare_y(prof::Profile, gene2num::Dict, numGene, numCell)
    fit_alpha = prof.var.fit_alpha
    train_y = Array{Float32,2}(undef, numGene, numCell)
    for (i, g) = enumerate(prof.var.index), j = 1:numCell
        k = get(gene2num, uppercase(g), 0)
        if k != 0
            train_y[k, j] = ismissing(fit_alpha[i]) ? zero(Float32) : fit_alpha[i]
        end
    end
    train_y
end

# A = adjacency_list(dg)
# train_X = prepare_X(prof, gene2num, d, numGene, numCell)
# train_y = prepare_y(prof, gene2num, numGene, numCell)
A = truncat_adjl!(adjacency_list(dg), numGene)
train_X = prepare_X(prof, truncat_gene2num!(gene2num, numGene), d, numGene, numCell)
train_y = prepare_y(prof, truncat_gene2num!(gene2num, numGene), numGene, numCell)

## Model

model = Chain(Concentration(d=>numGene),
              GeneRegulatory(A, numGene))

## Loss

loss(X, y) = mse(model(X), y)

function test_model(X, Y)
    Ŷ = model(X)
    @show size(Ŷ)
    @show loss(X, Y)
end

function test_gradient(X, Y)
    @show gradient(x -> sum(model(x)), X)
    @show gradient(x -> loss(x, Y), X)
end

