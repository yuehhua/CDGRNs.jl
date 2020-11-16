function load_data(dir::String; expr_file="gene_expr.tsv", var_file="genes.tsv", obs_file="obs.tsv")
    X = CSV.File(joinpath(dir, expr_file)) |> DataFrame
    var = CSV.File(joinpath(dir, var_file)) |> DataFrame
    obs = CSV.File(joinpath(dir, obs_file)) |> DataFrame
    Profile(SparseMatrixCSC(Matrix(X[:,2:end])'), obs, var)
end

function add_unspliced_data!(prof::Profile, dir::String;
                             unspliced_file="unspliced.tsv", spliced_file="spliced.tsv")
    unspliced = CSV.File(joinpath(dir, unspliced_file)) |> DataFrame
    spliced = CSV.File(joinpath(dir, spliced_file)) |> DataFrame
    prof.layers[:unspliced] = Matrix(unspliced[:,2:end])'
    prof.layers[:spliced] = Matrix(spliced[:,2:end])'
    prof
end

function add_velocity!(prof::Profile, dir::String;
    velo_file="velocity.tsv", var_velo_file="variance_velocity.tsv", velo_u_file="velocity_u.tsv")
    velocity = CSV.File(joinpath(dir, velo_file)) |> DataFrame
    variance_velocity = CSV.File(joinpath(dir, var_velo_file)) |> DataFrame
    velocity_u = CSV.File(joinpath(dir, velo_u_file)) |> DataFrame
    prof.layers[:velocity] = Matrix(velocity[:,2:end])'
    prof.layers[:var_velocity] = Matrix(variance_velocity[:,2:end])'
    prof.layers[:velocity_u] = collect(Missings.replace(Matrix(velocity_u[:,2:end])', 0))
    prof
end