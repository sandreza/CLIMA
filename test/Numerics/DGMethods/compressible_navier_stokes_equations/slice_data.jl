# rework data into slices
filepath = "baroclinic_slices.jld2"
gridhelper = GridHelper(cpu_grid.numerical)
ϕ = ScalarField(ρθ, gridhelper)


file = jldopen(filepath, "a+")
file["grid"] = cpu_grid
JLD2.Group(file, "state")
JLD2.Group(file, "time")
file["slices"][string(steps)] = Array(Q.ρθ)
file["time"][string(steps)] = time

