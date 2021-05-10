using MLDataUtils

function validate_score(reg::Type{MixtureRegression{K}}, X::AbstractVecOrMat, y::AbstractVector{T};
                        cv=0, λ=2e-2, criterion=aic, return_model=false)  where {K,T}
    if cv == 0
        model = fit(reg, X, y)
        score = criterion(model, λ=λ)
        return return_model ? Dict(:model=>model, :score=>score) : score
    else
        scores = T[]
        models = reg[]
        for (trains, vals) in kfolds(collect(1:length(y)), k=cv)
            model = fit(reg, _view(X, trains), view(y, trains))
            score = criterion(model, _view(X, vals), view(y, vals), λ=λ)
            return_model && push!(models, model)
            push!(scores, score)
        end
        return return_model ? Dict(:model=>models, :score=>scores) : scores
    end
end

select_hyperparams(scores, ::typeof(aic)) = argmin(scores)
select_hyperparams(scores, ::typeof(bic)) = argmin(scores)
select_hyperparams(scores, ::typeof(likelihood)) = argmax(scores)
select_hyperparams(scores, ::typeof(loglikelihood)) = argmax(scores)
select_hyperparams(scores, ::typeof(nll)) = argmin(scores)

function grid_search(reg::Type{MixtureRegression}, X::AbstractVecOrMat, y::AbstractVector{T}, k_range;
                     cv=0, λ=2e-2, criterion=aic, best_model=false, return_score=false, verbosity::Integer=0) where {T}
    mean_scores = T[]
    for k = k_range
        scores = validate_score(MixtureRegression{k}, X, y; cv=cv, λ=λ, return_model=false)
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
