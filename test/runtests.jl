using CDGRNs
using DataFrames
using CSV
using SnowyOwl
using Test

tests = [
    "data",
    "regression",
    "metrics",
    # "mixture",
]

include("test_utils.jl")

@testset "GRN" begin
    for t in tests
        include("$(t).jl")
    end
end
