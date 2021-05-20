function gmm_clustering(k::Integer, train::AbstractMatrix)
	gmm = GaussianMixtures.GMM(k, train, kind=:full)
    @warn "GMM(k=$k, n=$(size(train, 1))) failed, fallback to k-means."
	ll = GaussianMixtures.llpg(gmm, train)
	clusters = vec(map(x -> x[2], argmax(ll, dims=2)))
    return clusters
end

function kmeans_clustering(k::Integer, train::AbstractMatrix)
    n = size(train, 1)
    n < k && (k = n)
    @warn "k-means(k=$k, n=$n)"
	res = kmeans(train', k; maxiter=200, display=:none)
	clusters = assignments(res)
    return clusters
end

clustering(k::Integer, xs::AbstractVector, y::AbstractVector; method=gmm_clustering) = clustering(method, k, hcat(xs, y))
clustering(k::Integer, X::AbstractMatrix, y::AbstractVector; method=gmm_clustering) = clustering(method, k, hcat(X', y))
clustering(method, k::Integer, train::AbstractMatrix) = method(k, train)

function clustering(method::typeof(gmm_clustering), k::Integer, train::AbstractMatrix)
    try
        return method(k, train)
    catch e
        @warn "GMM(k=$k, n=$(size(train, 1))) failed, fallback to k-means."
    finally
        return kmeans_clustering(k, train)
    end
end
