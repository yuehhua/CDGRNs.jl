@model DPGLM(X, y, K) = begin
    D, N = size(X)
    
    α = 1.0
    w ~ Dirichlet(K, α)

	σ = Vector(undef, K)
	β₀ = Vector(undef, K)
	β = Matrix(undef, D, K)
	for i in 1:K
		σ[i] ~ InverseGamma(2, 3)
		β₀[i] ~ Normal(0, sqrt(3))
		β[:,i] ~ MvNormal(D, sqrt(10))
	end
    
    k = Vector{Int}(undef, N)
    for i in 1:N
        k[i] ~ Categorical(w)
    end
		
	σs = [σ[k[i]] for i in 1:N]
	μs = [β₀[k[i]] + X[:,i]' * β[:,k[i]] for i in 1:N]
	
	y ~ MvNormal(μs, sqrt.(σs))
    return k
end

# x isa Vector
@model DPGLM(x, y, K) = begin
    N = size(x, 1)
    
    α = 1.0
    w ~ Dirichlet(K, α)

	σ = Vector(undef, K)
	β₀ = Vector(undef, K)
	β = Vector(undef, K)
	for i in 1:K
		σ[i] ~ InverseGamma(2, 3)
		β₀[i] ~ Normal(0, sqrt(3))
		β[i] ~ Normal(0, sqrt(10))
	end
    
    k = Vector{Int}(undef, N)
    for i in 1:N
        k[i] ~ Categorical(w)
    end
		
	σs = [σ[k[i]] for i in 1:N]
	μs = [β₀[k[i]] + x[i] * β[k[i]] for i in 1:N]
	
	y ~ MvNormal(μs, sqrt.(σs))
    return k
end