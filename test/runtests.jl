using GRN
using GeometricFlux
using Flux
using GraphSignals
using DataFrames
using CSV
using SnowyOwl
using Zygote
using Test

tests = [
    "layers",
    "data",
    "regression",
    "metrics",
    "mixture",
]

include("test_utils.jl")

@testset "GRN" begin
    for t in tests
        include("$(t).jl")
    end
end
