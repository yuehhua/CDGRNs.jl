dir = joinpath(@__DIR__, "..", "results")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)


