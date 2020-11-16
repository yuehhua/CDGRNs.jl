using GRN
using GeometricFlux
using Flux
using GraphSignals
using DataFrames
using CSV
using SnowyOwl
using Test

tests = [
    "grn",
    "data"
]

@testset "GRN" begin
    for t in tests
        include("$(t).jl")
    end
end
