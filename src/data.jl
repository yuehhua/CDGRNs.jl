function load_data(dir::String; expr_file="gene_expr.tsv", var_file="genes.tsv", obs_file="obs.tsv")
    X = CSV.read(joinpath(dir, expr_file), DataFrame)
    var = CSV.read(joinpath(dir, var_file), DataFrame)
    obs = CSV.read(joinpath(dir, obs_file), DataFrame)
    Profile(SparseMatrixCSC(Matrix(X[:,2:end])'), var, obs)
end

function add_unspliced_data!(prof::Profile, dir::String;
                             unspliced_file="unspliced.tsv", spliced_file="spliced.tsv")
    unspliced = CSV.read(joinpath(dir, unspliced_file), DataFrame)
    spliced = CSV.read(joinpath(dir, spliced_file), DataFrame)
    prof.layers[:unspliced] = collect(Matrix(unspliced[:,2:end])')
    prof.layers[:spliced] = collect(Matrix(spliced[:,2:end])')
    prof
end

function add_moments!(prof::Profile, dir::String; unspliced_file="Mu.tsv", spliced_file="Ms.tsv")
    Mu = CSV.read(joinpath(dir, unspliced_file), DataFrame)
    Ms = CSV.read(joinpath(dir, spliced_file), DataFrame)
    prof.layers[:Mu] = collect(Matrix(Mu[:,2:end])')
    prof.layers[:Ms] = collect(Matrix(Ms[:,2:end])')
    prof
end

function add_velocity!(prof::Profile, dir::String;
    velo_file="velocity.tsv", var_velo_file="variance_velocity.tsv", velo_u_file="velocity_u.tsv")
    velocity = CSV.read(joinpath(dir, velo_file), DataFrame)
    variance_velocity = CSV.read(joinpath(dir, var_velo_file), DataFrame)
    velocity_u = CSV.read(joinpath(dir, velo_u_file), DataFrame)
    prof.layers[:velocity] = collect(Matrix(velocity[:,2:end])')
    prof.layers[:var_velocity] = collect(Matrix(variance_velocity[:,2:end])')
    prof.layers[:velocity_u] = collect(Missings.replace(Matrix(velocity_u[:,2:end])', 0))
    prof
end

function select_high_likelihood!(prof::Profile; min_likelihood=0.1)
    select_likelihood = x -> !ismissing(x) && x .≥ min_likelihood
    return filter!(:fit_likelihood => select_likelihood, prof)
end

function select_genes!(prof::Profile, gene_set)
    isin_geneset = x -> uppercase(x) in gene_set
    return filter!(:index => isin_geneset, prof)
end

load_tfs(filepath::String) = load(filepath, "tf_set")

function load_CHEA(dirpath::String)
    data_path = joinpath(dirpath, "data")
    filepath = joinpath(data_path, "gene_attribute_edges.txt")
    data = CSV.File(filepath, delim="\t") |> DataFrame
    regulations = DataFrame()
    regulations.tf = data.target[2:end]
    regulations.target = data.source[2:end]
    return regulations
end
