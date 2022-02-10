using CSV
using DataFrames

working_path = @__DIR__
data_path = joinpath(working_path, "data")

read_regulations(filepath) = CSV.File(filepath, delim="\t") |> DataFrame

function read_regulatory_pairs()
    data = read_regulations(joinpath(data_path, "gene_attribute_edges.txt"))
    regulations = DataFrame()
    regulations.tf = data.target[2:end]
    regulations.target = data.source[2:end]
    return regulations
end

regulations = read_regulatory_pairs()

tf_set = Set(regulations.tf)  # 199
target_set = Set(regulations.target)  # 21585

tf_set ∩ target_set  # 198
