module Plots

using Colors
using DataFrames
using MultivariateStats
using Plots, StatsPlots
using PyCall
gr()

const PROJECT_PATH = dirname(dirname(@__DIR__))

export 
    # tests
    test_pmf,
    test_cdf

include("utils.jl")
include("tree.jl")
include("regulations.jl")
include("visualization.jl")
include("tests.jl")
    
end
