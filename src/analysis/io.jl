function save_cdgrn(filename, cdgrn, tfs, prof, cntx)
    context_cor = cor(cdgrn, tfs, prof, cntx)
    context_cor = context_cor[.!isnan.(context_cor.ρ), :]
    context_cor.reg_type = context_cor.ρ .> 0
    context_cor.reg_stng = abs.(context_cor.ρ)
    CSV.write(filename, context_cor)
    return context_cor
end

function save_effective_gene_set(filename, context_cor, reg_strength)
    correlated = context_cor[context_cor.reg_stng .≥ reg_strength, :]
    gene_set = DataFrame(gene=unique!(vcat(correlated.tf, correlated.target)))
    CSV.write(filename, gene_set)
end
