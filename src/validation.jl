using MLDataUtils

function validate_score(reg::Type{MixtureRegression{K}}, X::AbstractVecOrMat, y::AbstractVector{T};
                        cv=0, λ=2e-2, criterion=aic, return_model=false)  where {K,T}
    if cv == 0
        init_cluster = clustering(K, X, y; method=gmm_clustering)
        model = fit(MixtureRegression, init_cluster, X, y)
        score = criterion(model, λ=λ)
        return return_model ? Dict(:model=>model, :score=>score) : score
    else
        scores = T[]
        models = MixtureRegression[]
        for (trains, vals) in kfolds(collect(1:length(y)), k=cv)
            init_cluster = clustering(K, _view(X, trains), view(y, trains); method=gmm_clustering)
            model = fit(MixtureRegression, init_cluster, _view(X, trains), view(y, trains))
            score = criterion(model, _view(X, vals), view(y, vals), λ=λ)
            return_model && push!(models, model)
            push!(scores, score)
        end
        return return_model ? Dict(:model=>models, :score=>scores) : scores
    end
end

function validate_score(reg::Type{GMR{K}}, X::AbstractVecOrMat{T};
                        λ=2e-2, criterion=aic, return_model=false)  where {K,T}
    model = fit(reg, X)
    score = criterion(model, X, λ=λ)
    return return_model ? (model=model, score=score, k=K) : (score=score, k=K)
end

select_hyperparams(scores, ::typeof(aic)) = argmin(scores)
select_hyperparams(scores, ::typeof(bic)) = argmin(scores)
select_hyperparams(scores, ::typeof(likelihood)) = argmax(scores)
select_hyperparams(scores, ::typeof(loglikelihood)) = argmax(scores)
select_hyperparams(scores, ::typeof(nll)) = argmin(scores)

lowest_score(::typeof(aic)) = Inf
lowest_score(::typeof(bic)) = Inf
lowest_score(::typeof(likelihood)) = -Inf
lowest_score(::typeof(loglikelihood)) = -Inf
lowest_score(::typeof(nll)) = -Inf

function grid_search(reg::Type{MixtureRegression}, X::AbstractVecOrMat, y::AbstractVector{T}, k_range;
                     cv=0, λ=2e-2, criterion=aic, best_model=false, return_score=false, verbosity::Integer=0) where {T}
    mean_scores = T[]
    for k = k_range
        scores = validate_score(reg{k}, X, y; cv=cv, λ=λ, return_model=false)
        ms = mean(scores)
        verbosity > 1 && @info "With k = $k"
        verbosity > 1 && @info "CV scores: $ms ($scores)"
        push!(mean_scores, ms)
    end

    best_k = k_range[select_hyperparams(mean_scores, criterion)]
    results = best_k

    if return_score
        results = Dict(:best_k=>best_k, :score=>mean_scores)
    end

    if best_model
        model = fit(MixtureRegression{best_k}, X, y)
        if return_score
            results[:model] = model
        else
            results = Dict(:best_k=>best_k, :model=>model)
        end
    end
    return results
end

function grid_search(reg::Type{GMR}, X::AbstractVecOrMat, k_range;
                     λ=2e-2, criterion=aic, verbosity::Integer=0)
    results = []
	for k = k_range
        res = validate_score(reg{k}, X; λ=λ, return_model=true)
        push!(results, res)
        verbosity > 1 && @info "With k = $(res[:k])"
        verbosity > 1 && @info "Score: $(res[:score])"
	end
    results
end

function best_result(results; criterion=aic)
    scores = map(x -> x[:score], results)
    return results[select_hyperparams(scores, criterion)]
end

function best_score(results; criterion=aic)
    scores = map(x -> x[:score], results)
    return scores[select_hyperparams(scores, criterion)]
end

function best_model(results; criterion=aic)
    scores = map(x -> x[:score], results)
    models = map(x -> x[:model], results)
    return models[select_hyperparams(scores, criterion)]
end
