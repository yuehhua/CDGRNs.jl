module Preprocess

using SparseArrays

using CSV
using DataFrames
using FileIO: load
using GaussianMixtureRegressions
using SnowyOwl
using OmicsProfiles

const PROJECT_PATH = dirname(dirname(@__DIR__))

export
    # data
    load_data,
    add_unspliced_data!,
    add_velocity!,
    add_moments!,
    load_tfs,
    load_CHEA,
    select_high_likelihood!,
    select_genes!,

    # pipeline
    load_profile,
    regulation_correlation,
    remove_spurious_pairs,

    # preprocess
    get_regulation_expr

include("data.jl")
include("pipeline.jl")
include("preprocess.jl")

end
