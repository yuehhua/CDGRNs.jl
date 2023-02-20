function test_pmf(ρ1, ρ2, condition1, condition2, plot_dir::String=@__DIR__; α=0.7, bincount=30,
        title="", save=[:svg, :png], figsize=(800, 600))
    df = DataFrame(
        ρ=vcat(ρ1, ρ2),
        condition=vcat(repeat([condition1], length(ρ1)), repeat([condition2], length(ρ2)))
    )
    min_ρ, max_ρ = extrema(df[!, :ρ])
    bins = range(min_ρ, stop=max_ρ, length=bincount)

    p = @df df histogram(:ρ, group=:condition, bins=bins,
        fillalpha=α, bar_edges=false,
        xlabel="TF-target pair correlation", ylabel="Count",
        thickness_scaling=2, widen=false, size=figsize,
    )

    SnowyOwl.Plots.save(p, joinpath(plot_dir, "histogram_$(title)"), save)
end

function test_cdf(ρ1, ρ2, condition1, condition2, plot_dir::String=@__DIR__; step=0.1,
        title="", save=[:svg, :png], figsize=(800, 600))
    r = -1.0:step:1.0
    cntx_cdf = DataFrame(x=r, y=ecdf(ρ1)(r), condition=condition1)
    global_cdf = DataFrame(x=r, y=ecdf(ρ2)(r), condition=condition2)
    df = vcat(cntx_cdf, global_cdf)

    p = @df df Plots.plot(:x, :y, group=:condition, linetype=:steppre,
        xlabel="TF-target pair correlation", ylabel="Cumulative count",
        yticks=[0., 0.25, 0.5, 0.75, 1.], legend_position=:bottomright,
        thickness_scaling=2, widen=false, size=figsize,
    )

    SnowyOwl.Plots.save(p, joinpath(plot_dir, "cdf_$(title)"), save)
end
