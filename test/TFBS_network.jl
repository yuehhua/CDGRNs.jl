using CSV, DataFrames
using JLD2
using GRN


# Load data
dir = joinpath("data", "fantom5_cat")
gene_tf = CSV.File(joinpath(dir, "gene_tf_map.csv")) |> DataFrame
first(gene_tf, 10)

tf = CSV.File(joinpath(dir, "tf_list.csv")) |> DataFrame
first(tf, 10)

tf_set = Set(tf.hgnc_id)  # 1990


# Make HGNC mapping to gene id
hgnc2id = make_mapping(tf, :hgnc_id=>:gene_id)

## not all tf has hgnc id
gene_tf = gene_tf[map(x -> x in tf_set, gene_tf.tf_hgnc_id), :]

gene_tf.tf_id = map(x -> hgnc2id[x], gene_tf.tf_hgnc_id)
first(gene_tf, 10)

gene_set = Set(gene_tf.gene_id)  # 42140
tf_set = Set(gene_tf.tf_id)  # 117


# Check gene set includes all tf set
issubset(tf_set, gene_set)  # false

gene_set = gene_set ∪ setdiff(tf_set, gene_set)


# Make sorted gene list
gene_set = sort([x for x = gene_set])
gene2num = Dict(x => i for (i, x) in enumerate(gene_set))


dg = make_graph(gene_tf, gene2num, length(gene_set))

@save "../results/tf_gene_network.jld2" dg
@save "../results/gene_set.jld2" gene_set