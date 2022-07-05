function context_correlation(tfs, prof, true_regulations, context, k)
    pairs = unique(zip(true_regulations.tf, true_regulations.target))
    context_cor = cor(tfs, prof, pairs, context)
    context_pairs = collect(zip(context_cor.tf, context_cor.target))
    context_cor[!, :context], selected_context = max_cor(context_cor, k)
    return context_cor, selected_context, context_pairs
end

test_pmf(ρ1, ρ2) = MannWhitneyUTest(abs.(ρ1), abs.(ρ2))

test_cdf(ρ1, ρ2) = ApproximateTwoSampleKSTest(ρ1, ρ2)
