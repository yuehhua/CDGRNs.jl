savefig(plt, path::String, ::Nothing) = nothing
savefig(plt, path::String, ext) = Plots.savefig(plt, "$(path).$(ext)")

function savefig(plt, path::String, exts::AbstractVector)
    for ext in exts
        Plots.savefig(plt, "$(path).$(ext)")
    end
end
