dir = joinpath(@__DIR__, "..", "results", "pancreas")
prof = load_data(dir)
add_unspliced_data!(prof, dir)
add_velocity!(prof, dir)


