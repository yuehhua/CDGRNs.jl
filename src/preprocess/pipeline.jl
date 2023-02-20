function load_profile(dir::String)
    prof = load_data(dir)
    add_unspliced_data!(prof, dir)
    add_velocity!(prof, dir)
    add_moments!(prof, dir)
    return prof
end

function regulation_correlation(filename)
    total_results = load(filename, "total_results")
    nonsingle_pairs = filter(x -> x[:best_k] != 1, total_results)
    cor_pairs = corr_table(nonsingle_pairs)
    cor_pairs.is_tf = cor_pairs.tf .== cor_pairs.target
    return cor_pairs, nonsingle_pairs
end

function remove_spurious_pairs(cor_pairs, nonsingle_pairs)
    # map to database
    database = load_CHEA(joinpath(PROJECT_PATH, "CHEA"))
    pairset = make_pairset(database)
    cor_pairs.is_regulation = query_pairset(cor_pairs, pairset)
    true_regulations = cor_pairs[cor_pairs.is_regulation .& .!cor_pairs.is_tf, :]
    true_reg_pairs = filter(x -> (uppercase(x[:tf_name]), uppercase(x[:gene_name])) in pairset, nonsingle_pairs)
    return true_regulations, true_reg_pairs
end
