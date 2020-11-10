using GRN
using GeometricFlux
using Flux
using GraphSignals
using Test

tests = [
    "grn"
]

@testset "GRN" begin
    for t in tests
        include("$(t).jl")
    end
end
