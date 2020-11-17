function make_mapping(df::DataFrame, p::Pair{Symbol,Symbol})
    Dict(df[i, p.first] => df[i, p.second] for i = 1:nrow(df))
end