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
    prof.layers[:unspliced] = Matrix(unspliced[:,2:end])'
    prof.layers[:spliced] = Matrix(spliced[:,2:end])'
    prof
end

function add_moments!(prof::Profile, dir::String; unspliced_file="Mu.tsv", spliced_file="Ms.tsv")
    Mu = CSV.read(joinpath(dir, unspliced_file), DataFrame)
    Ms = CSV.read(joinpath(dir, spliced_file), DataFrame)
    prof.layers[:Mu] = Matrix(Mu[:,2:end])'
    prof.layers[:Ms] = Matrix(Ms[:,2:end])'
    prof
end

function add_velocity!(prof::Profile, dir::String;
    velo_file="velocity.tsv", var_velo_file="variance_velocity.tsv", velo_u_file="velocity_u.tsv")
    velocity = CSV.read(joinpath(dir, velo_file), DataFrame)
    variance_velocity = CSV.read(joinpath(dir, var_velo_file), DataFrame)
    velocity_u = CSV.read(joinpath(dir, velo_u_file), DataFrame)
    prof.layers[:velocity] = Matrix(velocity[:,2:end])'
    prof.layers[:var_velocity] = Matrix(variance_velocity[:,2:end])'
    prof.layers[:velocity_u] = collect(Missings.replace(Matrix(velocity_u[:,2:end])', 0))
    prof
end

# function filter_genes!(prof::Profile; min_likelihood=0.1)
#     select_likelihood = x -> !ismissing(x) && x .≥ min_likelihood

#     vars = filter(:fit_likelihood => select_likelihood, prof.var)
#     prof.var = sort!(vars, :fit_likelihood, rev=true)
    
#     selected_rows = select_likelihood.(prof.var.fit_likelihood)
#     prof.data = prof.data[selected_rows, :]
#     prof.layers[:Mu] = prof.layers[:Mu][selected_rows, :]
#     prof.layers[:velocity_u] = prof.layers[:velocity_u][selected_rows, :]
#     prof.layers[:Ms] = prof.layers[:Ms][selected_rows, :]
#     prof.layers[:velocity] = prof.layers[:velocity][selected_rows, :]
#     prof
# end

load_tfs(filepath::String) = load(filepath, "tf_set")

# function filter_tfs!(prof::Profile, tf_set)
#     select_index = x -> uppercase(x) in tf_set
#     select_likelihood = x -> !ismissing(x)

#     selected_rows = select_index.(prof.var.index)
#     tf_vars = filter(:index => select_index, prof.var)
#     tf_data = prof.data[selected_rows, :]
#     tf_u = prof.layers[:Mu][selected_rows, :]
#     tf_vᵤ = prof.layers[:velocity_u][selected_rows, :]
#     tf_s = prof.layers[:Ms][selected_rows, :]
#     tf_vₛ = prof.layers[:velocity][selected_rows, :]

#     selected_rows = select_likelihood.(tf_vars.fit_likelihood)
#     filter!(:fit_likelihood => x -> select_likelihood(x), tf_vars)
#     tf_data = tf_data[selected_rows, :]
#     tf_u = tf_u[selected_rows, :]
#     tf_vᵤ = tf_vᵤ[selected_rows, :]
#     tf_s = tf_s[selected_rows, :]
#     tf_vₛ = tf_vₛ[selected_rows, :]
# end
