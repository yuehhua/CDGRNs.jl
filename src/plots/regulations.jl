function regulations(data::DataFrame, x::String, y::String, z::String; model=nothing, spliced=false,
                     figsize=(1000, 600))
    x_sym = Symbol(x * "_s")
    y_sym = Symbol(y * "_s")
    z_sym = spliced ? Symbol(z * "_s") : Symbol(z * "_u")
    xlabel = "Spliced mRNA of $x"
    ylabel = "Spliced mRNA of $y"
    zlabel = spliced ? "Spliced mRNA of $z" : "Unspliced mRNA of $z"
    p = scatter(data[:,x_sym], data[:,y_sym], data[:,z_sym], group=data[:,:cell_type],
                xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                markersize=2, markerstrokewidth=0, legend_position=:outerright,
                thickness_scaling=2, widen=false, dpi=300, size=figsize,
        )
    if !isnothing(model)
        x = range(floor(minimum(data[:,x_sym])), stop=ceil(maximum(data[:,x_sym])), length=100)
        y = range(floor(minimum(data[:,y_sym])), stop=ceil(maximum(data[:,y_sym])), length=100)
        plot!(p, x, y, (x,y) -> predict(model, [x, y]'), st=:surface,
            c=:grey, fillalpha=0.7, cb=false)
    end
    return p
end

function regulations(data::DataFrame, x::String, y::String; model=nothing, spliced=false,
                     figsize=(1000, 600))
    x_sym = Symbol(x * "_s")
    y_sym = spliced ? Symbol(y * "_s") : Symbol(y * "_u")
    xlabel = "Spliced mRNA of $x"
    ylabel = spliced ? "Spliced mRNA of $y" : "Unspliced mRNA of $y"
    p = scatter(data[:,x_sym], data[:,y_sym], group=data[:,:cell_type],
                xlabel=xlabel, ylabel=ylabel,
                markersize=3, markerstrokewidth=0, legend_position=:outerright,
                thickness_scaling=2, widen=false, dpi=300, size=figsize,
        )
    if !isnothing(model)
        x = range(floor(minimum(data[:,x_sym])), stop=ceil(maximum(data[:,x_sym])), length=100)
        plot!(p, x, y -> predict(model, [y]'), st=:surface,
            c=:grey, fillalpha=0.7, cb=false)
    end
    return p
end
