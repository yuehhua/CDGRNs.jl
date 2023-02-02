module Plots

using Colors
using DataFrames
using GaussianMixtureRegressions
using MultivariateStats
using Plots, StatsPlots
using PyCall
using UMAP
gr()

const PROJECT_PATH = dirname(dirname(@__DIR__))

include("save.jl")

include("utils.jl")
include("tree.jl")
include("regulations.jl")
include("visualization.jl")
include("tests.jl")
include("landscapes.jl")

end
