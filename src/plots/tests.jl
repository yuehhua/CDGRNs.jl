function test_pmf(ρ1, ρ2, condition1, condition2, plot_dir::String=""; α=0.7, bincount=30,
        title="", save_png=true, save_svg=true, figsize=(800, 600))
    ρ = vcat(ρ1, ρ2)
    condition = vcat(repeat([condition1], length(ρ1)), repeat([condition2], length(ρ2)))
    plot_df = DataFrame(ρ=ρ, condition=condition)
    bins = range(minimum(plot_df[!, :ρ]), stop=maximum(plot_df[!, :ρ]), length=bincount)

    p = @df plot_df histogram(:ρ, group=:condition, bins=bins,
        fillalpha=α, bar_edges=false,
        xlabel="TF-target pair correlation", ylabel="Count",
        thickness_scaling=2, widen=false, size=figsize,
    )
    save_svg && savefig(p, joinpath(plot_dir, "histogram_$(title).svg"))
    save_png && savefig(p, joinpath(plot_dir, "histogram_$(title).png"))
end

function test_cdf(ρ1, ρ2, condition1, condition2, plot_dir::String=""; step=0.1,
        title="", save_png=true, save_svg=true, figsize=(800, 600))
    r = -1.0:step:1.0
    cntx_cdf = DataFrame(x=r, y=ecdf(ρ1)(r), condition=condition1)
    global_cdf = DataFrame(x=r, y=ecdf(ρ2)(r), condition=condition2)
    cdf_df = vcat(cntx_cdf, global_cdf)

    p = @df cdf_df Plots.plot(:x, :y, group=:condition, linetype=:steppre,
        xlabel="TF-target pair correlation", ylabel="Cumulative count",
        yticks=[0., 0.25, 0.5, 0.75, 1.], legend_position=:bottomright,
        thickness_scaling=2, widen=false, size=figsize,
    )
    save_svg && savefig(p, joinpath(plot_dir, "cdf_$(title).svg"))
    save_png && savefig(p, joinpath(plot_dir, "cdf_$(title).png"))
end
