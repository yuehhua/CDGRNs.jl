using CDGRNs
using CSV, DataFrames
using SnowyOwl
using Test

tests = [
    # "data",
]

@testset "GRN" begin
    for t in tests
        include("$(t).jl")
    end
end
