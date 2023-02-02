function load_data(dir::String; expr_file="gene_expr.tsv", var_file="genes.tsv", obs_file="obs.tsv")
    X = CSV.read(joinpath(dir, expr_file), DataFrame)
    var = CSV.read(joinpath(dir, var_file), DataFrame)
    obs = CSV.read(joinpath(dir, obs_file), DataFrame)
    return Profile(SparseMatrixCSC{Float64, UInt32}(Matrix(X[:,2:end])'), :RNA, var, obs,
        varindex=:index, obsindex=:index)
end

function add_unspliced_data!(prof::AnnotatedProfile, dir::String;
                             unspliced_file="unspliced.tsv", spliced_file="spliced.tsv")
    unspliced = CSV.read(joinpath(dir, unspliced_file), DataFrame)
    spliced = CSV.read(joinpath(dir, spliced_file), DataFrame)
    op = prof.omics[:RNA]
    op.layers[:unspliced] = collect(Matrix(unspliced[:,2:end])')
    op.layers[:spliced] = collect(Matrix(spliced[:,2:end])')
    return prof
end

function add_moments!(prof::AnnotatedProfile, dir::String; unspliced_file="Mu.tsv", spliced_file="Ms.tsv")
    Mu = CSV.read(joinpath(dir, unspliced_file), DataFrame)
    Ms = CSV.read(joinpath(dir, spliced_file), DataFrame)
    op = prof.omics[:RNA]
    op.layers[:Mu] = collect(Matrix(Mu[:,2:end])')
    op.layers[:Ms] = collect(Matrix(Ms[:,2:end])')
    prof
end

function add_velocity!(prof::AnnotatedProfile, dir::String;
    velo_file="velocity.tsv", var_velo_file="variance_velocity.tsv", velo_u_file="velocity_u.tsv")
    velocity = CSV.read(joinpath(dir, velo_file), DataFrame)
    variance_velocity = CSV.read(joinpath(dir, var_velo_file), DataFrame)
    velocity_u = CSV.read(joinpath(dir, velo_u_file), DataFrame)
    op = prof.omics[:RNA]
    op.layers[:velocity] = collect(Matrix(velocity[:,2:end])')
    op.layers[:var_velocity] = collect(Matrix(variance_velocity[:,2:end])')
    op.layers[:velocity_u] = collect(Missings.replace(Matrix(velocity_u[:,2:end])', 0))
    return prof
end

function select_high_likelihood!(prof::AnnotatedProfile; min_likelihood=0.1, omics=:RNA)
    select_likelihood = x -> !ismissing(x) && x .≥ min_likelihood
    filter!(:fit_likelihood => select_likelihood, prof.omics[omics])
    return prof
end

select_genes(prof::AnnotatedProfile, gene_set) = select_genes!(copy(prof), gene_set)

function select_genes!(prof::AnnotatedProfile, gene_set; omics=:RNA)
    isin_geneset = x -> uppercase(x) in gene_set
    filter!(:index => isin_geneset, prof.omics[omics])
    return prof
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

function get_regulation_expr(prof::AnnotatedProfile, tfs::AnnotatedProfile, true_reg::DataFrame;
            labels=nothing, latent_time=false)
    tf_list = unique(true_reg.tf)
    gene_list = unique(true_reg.target)

    df = DataFrame()
    isnothing(labels) || (df[!, :cell] = prof.obs[!, labels])
    latent_time && (df[!, :time] = prof.obs.latent_time)
    for gene in gene_list
        df[!, Symbol(gene * "_u")] = vec(prof.RNA.Mu[gene, :])
    end
    for tf in tf_list
        df[!, Symbol(tf * "_s")] = vec(tfs.RNA.Ms[tf, :])
    end
    return df
end
