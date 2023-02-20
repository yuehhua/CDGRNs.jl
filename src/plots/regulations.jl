function regulations(data::DataFrame, x::String, y::String, z::String; model=nothing, spliced=false,
                     figsize=(1000, 600))
    x_sym = Symbol(x * "_s")
    y_sym = Symbol(y * "_s")
    z_sym = spliced ? Symbol(z * "_s") : Symbol(z * "_u")
    xlabel = "Spliced mRNA of $x"
    ylabel = "Spliced mRNA of $y"
    zlabel = spliced ? "Spliced mRNA of $z" : "Unspliced mRNA of $z"
    xs, ys, zs = data[:,x_sym], data[:,y_sym], data[:,z_sym]
    p = scatter(xs, ys, zs; group=data[:,:cell_type],
        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
        markersize=2, figsize=figsize)
    if !isnothing(model)
        x = range(floor(minimum(xs)), stop=ceil(maximum(xs)), length=100)
        y = range(floor(minimum(ys)), stop=ceil(maximum(ys)), length=100)
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
    xs, ys = data[:,x_sym], data[:,y_sym]
    p = scatter(xs, ys; group=data[:,:cell_type],
        xlabel=xlabel, ylabel=ylabel, markersize=3, figsize=figsize)
    if !isnothing(model)
        x = range(floor(minimum(xs)), stop=ceil(maximum(xs)), length=100)
        plot!(p, x, y -> predict(model, [y]'), st=:surface,
            c=:grey, fillalpha=0.7, cb=false)
    end
    return p
end
