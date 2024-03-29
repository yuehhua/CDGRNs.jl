module Plots

using Colors
using DataFrames
using GaussianMixtureRegressions
using Plots, StatsPlots
using PyCall
gr()

const PROJECT_PATH = dirname(dirname(@__DIR__))

include("utils.jl")
include("tree.jl")
include("regulations.jl")
include("visualization.jl")
include("tests.jl")
include("landscapes.jl")

end
