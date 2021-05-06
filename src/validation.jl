using MLDataUtils

function cross_val_score(reg::Type{MixtureRegression{K}}, X::AbstractVecOrMat, y::AbstractVector{T};
                         cv=5, λ=2e-2, criterion=aic)  where {K,T}
    scores = T[]
    for (trains, vals) in kfolds(collect(1:length(y)), k=cv)
        model = fit(reg, _view(X, trains), view(y, trains))
        score = criterion(model, _view(X, vals), view(y, vals), λ=λ)
        push!(scores, score)
    end
    return scores
end

function grid_search(reg::Type{MixtureRegression}, X::AbstractVecOrMat, y::AbstractVector{T}, k_range;
                     cv=5, λ=2e-2, criterion=aic, verbosity::Integer=0) where {T}
    mean_scores = T[]
    for k = k_range
        scores = cross_val_score(MixtureRegression{k}, X, y; cv=cv, λ=λ)
        ms = mean(scores)
        verbosity > 1 && @info "With k = $k"
        verbosity > 1 && @info "CV scores: $ms ($scores)"
        push!(mean_scores, ms)
    end
    return k_range[argmin(mean_scores)]
end