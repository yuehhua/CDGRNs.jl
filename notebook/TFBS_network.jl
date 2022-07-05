using CSV, DataFrames
using JLD2, FileIO
using CDGRNs


# Load data
dir = joinpath("data", "fantom5_cat")
gene_tf = CSV.File(joinpath(dir, "gene_tf_map.csv")) |> DataFrame
first(gene_tf, 10)

tf = CSV.File(joinpath(dir, "tf_list.csv")) |> DataFrame
first(tf, 10)

hgnc = CSV.File(joinpath(dir, "hgnc_info.csv")) |> DataFrame

tf_set = Set(tf.hgnc_id)  # 1990

# To hgnc id
ensemble2hgnc = CDGRNs.make_mapping(hgnc, :ensembl_gene_id=>:hgnc_id)

gene_ids = map(x -> split(x, '.')[1], gene_tf.gene_id)
gene_tf.gene_hgnc = map(x -> get(ensemble2hgnc, x, missing), gene_ids)

gene_tf = dropmissing(gene_tf)


# Select by gene list
gene_list = CSV.File(joinpath(dir, "gene_list.csv")) |> DataFrame

hgnc = hgnc[map(x -> x in gene_list.symbol, hgnc.symbol), :]

# To gene symbol
ensemble2symb = CDGRNs.make_mapping(hgnc, :ensembl_gene_id=>:symbol)
hgnc2symb = CDGRNs.make_mapping(hgnc, :hgnc_id=>:symbol)

gene_ids = map(x -> split(x, '.')[1], gene_tf.gene_id)
gene_tf.gene_symbol = map(x -> get(ensemble2symb, x, missing), gene_ids)
gene_tf.tf_symbol = map(x -> get(hgnc2symb, x, missing), gene_tf.tf_hgnc_id)


# Select mapped entities only
gene_tf = dropmissing(gene_tf)
first(gene_tf, 10)

gene_set = Set(gene_tf.gene_symbol)  # 4569
tf_set = Set(gene_tf.tf_symbol)  # 36

# Union all tfs into gene set
gene_set = gene_set ∪ tf_set


# Make sorted gene list
gene_set = sort([x for x = gene_set])
tf_set = sort([x for x = tf_set])
gene2num = Dict(x => i for (i, x) in enumerate(gene_set))

# reverse arrow means genes are regulated by TFs
dg = CDGRNs.make_graph(gene_tf, gene2num, length(gene_set))

save("results/tf_gene_network.jld2", "dg", dg)
save("results/gene_set.jld2", "gene_set", gene_set)
save("results/tf_set.jld2", "tf_set", tf_set)
save("results/gene2num.jld2", "gene2num", gene2num)
