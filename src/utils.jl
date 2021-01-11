function make_mapping(df::DataFrame, p::Pair{Symbol,Symbol})
    Dict(df[i, p.first] => df[i, p.second] for i = 1:nrow(df))
end

function truncat_adjl!(adjl::Vector, top_n)
    adjl = adjl[1:top_n]
    for i = 1:top_n
        l = adjl[i]
        adjl[i] = l[l .<= top_n]
    end
    adjl
end

function truncat_gene2num!(gene2num::Dict, top_n)
    for (k, v) in gene2num
        if v > top_n
            delete!(gene2num, k)
        end
    end
    gene2num
end
