struct DistributionTransformation{T<:Distribution,S<:Distribution}
    src::T
    dst::S
end

function fit(::Type{DistributionTransformation}, dist::Type{<:Distributions.Distribution}, xs::AbstractArray)
	model = Distributions.fit(dist, xs)
	normal = Normal(mean(model), std(model))
    DistributionTransformation(model, normal)
end

function transform!(trans::DistributionTransformation, xs::AbstractArray)
    for i = CartesianIndices(xs)
        xs[i] = quantile(trans.dst, cdf(trans.src, xs[i]))
    end
end

function transform(trans::DistributionTransformation, xs::AbstractArray)
    quantile.(trans.dst, cdf.(trans.src, xs))
end

function fit_transform(::Type{DistributionTransformation}, dist::Type{<:Distributions.Distribution}, xs::AbstractArray)
    dt = fit(DistributionTransformation, dist, xs)
    transform(dt, xs)
end
