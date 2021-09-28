using CSV, DataFrames
using LightGraphs
using JLD2
using FileIO
using CDGRN
using LibPQ

const RESULT_PATH = joinpath(CDGRN.PROJECT_PATH, "results")
const DATA_PATH = joinpath(CDGRN.PROJECT_PATH, "data", "fantom5_cat")

# Create TFBS network

## Load data
gene_tf = CSV.File(joinpath(DATA_PATH, "gene_tf_map.csv")) |> DataFrame
tf = CSV.File(joinpath(DATA_PATH, "tf_list.csv")) |> DataFrame
total_tf_set = Set(tf.hgnc_id)

## Make HGNC mapping to gene id
hgnc2id = Dict(tf[i, :hgnc_id] => tf[i, :gene_id] for i = 1:nrow(tf))

# not all TF has hgnc id
filter!(:tf_hgnc_id => x -> x in total_tf_set, gene_tf)

gene_tf.tf_id = map(x -> hgnc2id[x], gene_tf.tf_hgnc_id)
hgnc_gene_set = Set(gene_tf.gene_id)
hgnc_tf_set = Set(gene_tf.tf_id)

## Check if gene set includes all tf set
hgnc_tf_set ⊆ hgnc_gene_set

## Make sorted gene list
hgnc_gene_tf_set = sort([x for x = hgnc_tf_set ∪ hgnc_gene_set])
gene2num = Dict(x => i for (i, x) in enumerate(hgnc_gene_tf_set))

## Make simple directed graph

dg = SimpleDiGraph(length(hgnc_gene_tf_set))
for i = 1:nrow(gene_tf)
	g = gene2num[gene_tf[i, :gene_id]]
	tf = gene2num[gene_tf[i, :tf_id]]
	add_edge!(dg, tf, g)
end

save(joinpath(RESULT_PATH, "tf_gene_network.jld2"), "dg", dg)
save(joinpath(RESULT_PATH, "gene_set.jld2"), "gene_set", hgnc_gene_tf_set)

## Read HGCN_INFO

function from_sql(sql)
	conn = LibPQ.Connection("dbname=fantom5_cat user=yuehhua host=localhost port=5432")
	result = execute(conn, sql)
	df = result |> DataFrame
	close(result)
	close(conn)
	return df
end

sql = "select hgnc_id, symbol, name from hgnc_info;"
hgnc_info = from_sql(sql)
gene_tf = leftjoin(gene_tf, hgnc_info, on=:tf_hgnc_id=>:hgnc_id, renamecols = "" => "_tf")
rename!(gene_tf, Dict(:symbol_tf=>:tf_symbol, :name_tf=>:tf_name))
tf_set = unique(skipmissing(gene_tf.tf_symbol))
save(joinpath(RESULT_PATH, "tf_set.jld2"), "tf_set", tf_set)
