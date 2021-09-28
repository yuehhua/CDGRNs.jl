### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 8d81f736-c99d-4f34-bfac-b1063b37f209
begin
	using Distributions
	using StatsPlots
	using Random
	using Turing
	using MCMCChains
	using RDatasets
	using MLDataUtils: shuffleobs, splitobs, rescale!
	using Distances
	Turing.setprogress!(false);
end

# ╔═╡ 5a5aaee8-1e77-463e-ba2f-82c00a94bf0d
html"""
<style>
	main {
		max-width: 90%;
	}
</style>
"""

# ╔═╡ a099df8c-d2ad-11eb-0f99-b151021a7a92
md"# Example of Gaussian mixture model"

# ╔═╡ 2a77ff15-ff3f-4e94-846f-0851617706b1
md"## Generate data"

# ╔═╡ 339eb4ca-9257-4353-bdac-7ef383fdad56
begin
	Random.seed!(3)
	N = 30

	# Parameters for each cluster, we assume that each cluster is Gaussian distributed in the example.
	μs = [-3.5, 0.0]
	x = mapreduce(c -> rand(MvNormal([μs[c], μs[c]], 1.), N), hcat, 1:2)
end

# ╔═╡ 1483029c-d686-41a9-82d7-44c22057f808
scatter(x[1,:], x[2,:], legend = false, title = "Synthetic Dataset")

# ╔═╡ 1c506a19-b439-451f-af0f-f032c0e3a4a9
md"## Gaussian Mixture Model in Turing"

# ╔═╡ 276af183-0bc4-4f51-b076-657408da52e0
@model GaussianMixtureModel(x, K) = begin
    D, N = size(x)
    
    α = 1.0
    w ~ Dirichlet(K, α)

	μ = Vector(undef, K)
	σ = Vector(undef, K)
	for i in 1:K
		σ[i] ~ InverseGamma(2,3)
		μ[i] ~ Normal(0, sqrt(σ[i]))
	end
    
    k = Vector{Int}(undef, N)
    for i in 1:N
        k[i] ~ Categorical(w)
		μs = [μ[k[i]] for j in 1:K]
		σs = [σ[k[i]] for j in 1:K]
        x[:,i] ~ MvNormal(μs, σs)
    end
    return k
end

# ╔═╡ 5a9a0b77-a1fd-46e1-b6a7-efd08f116b7a
begin
	K = 2
	gmm_model = GaussianMixtureModel(x, K);
end

# ╔═╡ 651743ff-7ed8-4863-9f39-2c8adae54572
begin
	gmm_sampler = Gibbs(PG(100, :k), HMC(0.05, 10, :μ))
	tchain = mapreduce(c -> sample(gmm_model, gmm_sampler, 100), chainscat, 1:3);
end

# ╔═╡ 82fdac4b-cd9a-4e3b-954c-549cb5c4b460
md"## Visualization"

# ╔═╡ 262b9ed9-d170-442f-8224-094cf5b2a21f
begin
	ids = findall(map(name -> occursin("μ", string(name)), names(tchain)));
	p=plot(tchain[:, ids, :], legend=true, labels = ["Mu 1" "Mu 2"], colordim=:parameter)
end

# ╔═╡ 14d0c3e4-41af-45a0-ac31-2d545505ed56
function predict(x, y, w, μ)
    # Use log-sum-exp trick for numeric stability.
    return Turing.logaddexp(
        log(w[1]) + logpdf(MvNormal([μ[1], μ[1]], 1.), [x, y]), 
        log(w[2]) + logpdf(MvNormal([μ[2], μ[2]], 1.), [x, y])
    )
end

# ╔═╡ 82a493d6-7c2e-4abd-8ffa-020f76aa136d
begin
	contour(range(-5, stop = 3), range(-6, stop = 2), 
		(x, y) -> predict(x, y, [0.5, 0.5], [mean(tchain[:, :, 3][Symbol("μ[1]")]), mean(tchain[:, :, 3][Symbol("μ[2]")])])
	)
	scatter!(x[1,:], x[2,:], legend = false, title = "Synthetic Dataset")
end

# ╔═╡ a82997d9-acff-42fa-b610-c2f3827028de
md"## Linear regression"

# ╔═╡ 6b473e62-1401-44a4-93de-c85cdc9164ea
md"## Load data"

# ╔═╡ a0686aba-ec1f-44e7-bafd-350409de95ba
begin
	data = RDatasets.dataset("datasets", "mtcars");
	
	# Remove the model column.
	select!(data, Not(:Model))

	# Split our dataset 70%/30% into training/test sets.
	trainset, testset = splitobs(shuffleobs(data), 0.7)

	# Turing requires data in matrix form.
	target = :MPG
	train = Matrix(select(trainset, Not(target)))
	test = Matrix(select(testset, Not(target)))
	train_target = trainset[:, target]
	test_target = testset[:, target]

	# Standardize the features.
	μ, σ = rescale!(train; obsdim = 1)
	rescale!(test, μ, σ; obsdim = 1)

	# Standardize the targets.
	μtarget, σtarget = rescale!(train_target; obsdim = 1)
	rescale!(test_target, μtarget, σtarget; obsdim = 1)
end

# ╔═╡ 9457ae6f-1193-414e-b537-8f3222b72948
# Bayesian linear regression.
@model function linear_regression(x, y)
    # Set variance prior.
    σ₂ ~ truncated(Normal(0, 100), 0, Inf)
    
    # Set intercept prior.
    intercept ~ Normal(0, sqrt(3))
    
    # Set the priors on our coefficients.
    nfeatures = size(x, 2)
    coefficients ~ MvNormal(nfeatures, sqrt(10))
    
    # Calculate all the mu terms.
    mu = intercept .+ x * coefficients
    y ~ MvNormal(mu, sqrt(σ₂))
end

# ╔═╡ d73df181-37a0-402e-b1e1-782a72e5307d
begin
	model = linear_regression(train, train_target)
	chain = sample(model, NUTS(0.65), 3_000);
end

# ╔═╡ 8324fe70-c254-460e-964b-522b51d208e8
md"## Visualize"

# ╔═╡ 6d5f3d46-674f-4d52-9ca5-414dbc12fe39
plot(chain)

# ╔═╡ b722ea69-f0f0-48b2-9046-55ca31b0b7a2
describe(chain)

# ╔═╡ Cell order:
# ╟─5a5aaee8-1e77-463e-ba2f-82c00a94bf0d
# ╟─a099df8c-d2ad-11eb-0f99-b151021a7a92
# ╠═8d81f736-c99d-4f34-bfac-b1063b37f209
# ╟─2a77ff15-ff3f-4e94-846f-0851617706b1
# ╠═339eb4ca-9257-4353-bdac-7ef383fdad56
# ╠═1483029c-d686-41a9-82d7-44c22057f808
# ╟─1c506a19-b439-451f-af0f-f032c0e3a4a9
# ╠═276af183-0bc4-4f51-b076-657408da52e0
# ╠═5a9a0b77-a1fd-46e1-b6a7-efd08f116b7a
# ╠═651743ff-7ed8-4863-9f39-2c8adae54572
# ╟─82fdac4b-cd9a-4e3b-954c-549cb5c4b460
# ╠═262b9ed9-d170-442f-8224-094cf5b2a21f
# ╠═14d0c3e4-41af-45a0-ac31-2d545505ed56
# ╠═82a493d6-7c2e-4abd-8ffa-020f76aa136d
# ╟─a82997d9-acff-42fa-b610-c2f3827028de
# ╟─6b473e62-1401-44a4-93de-c85cdc9164ea
# ╠═a0686aba-ec1f-44e7-bafd-350409de95ba
# ╠═9457ae6f-1193-414e-b537-8f3222b72948
# ╠═d73df181-37a0-402e-b1e1-782a72e5307d
# ╟─8324fe70-c254-460e-964b-522b51d208e8
# ╠═6d5f3d46-674f-4d52-9ca5-414dbc12fe39
# ╠═b722ea69-f0f0-48b2-9046-55ca31b0b7a2
