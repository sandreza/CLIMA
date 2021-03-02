# rework data into slices
filepath = "baroclinic_states.jld2"
readfile = jldopen(filepath, "r+")
Nx = 128
Ny = 256
Nz = 8
filepath = "slices.jld2"
gridhelper = GridHelper(cpu_grid.numerical)
ϕ = ScalarField(ρθ[:,1,:], gridhelper)

Ωˣ = cpu_grid.domain[1]
Ωʸ = cpu_grid.domain[2]
Ωᶻ = cpu_grid.domain[3]

file = jldopen(filepath, "a+")
file["grid"] = cpu_grid
JLD2.Group(file, "state")
JLD2.Group(file, "time")

for key in keys(readfile["state"])
    println(key)
    rρθ =  readfile["state"][key][:,1,:]
    file["time"][key] = readfile["time"][key]

    ϕ .= rρθ

    # edge 1
    x = range(Ωˣ.min, Ωˣ.max, length = Nx)
    y = range(Ωʸ.min, Ωʸ.min, length = 1)
    z = range(Ωᶻ.min, Ωᶻ.max, length = Nz)
    ϕedge1 = ϕ(x, y, z)[:,1,:]
    # edge 2
    x = range(Ωˣ.min, Ωˣ.max, length = Nx)
    y = range(Ωʸ.max, Ωʸ.max, length = 1)
    z = range(Ωᶻ.min, Ωᶻ.max, length = Nz)
    ϕedge2 = ϕ(x, y, z)[:,1,:]
    # edge 3
    x = range(Ωˣ.min, Ωˣ.min, length = 1)
    y = range(Ωʸ.min, Ωʸ.max, length = Ny)
    z = range(Ωᶻ.min, Ωᶻ.max, length = Nz)
    ϕedge3 = ϕ(x, y, z)[1,:,:]
    # edge 4
    x = range(Ωˣ.max, Ωˣ.max, length = 1)
    y = range(Ωʸ.min, Ωʸ.max, length = Ny)
    z = range(Ωᶻ.min, Ωᶻ.max, length = Nz)
    ϕedge4 = ϕ(x, y, z)[1,:,:]
    # edge 5
    x = range(Ωˣ.min, Ωˣ.max, length = Nx)
    y = range(Ωʸ.min, Ωʸ.max, length = Ny)
    z = range(Ωᶻ.min, Ωᶻ.min, length = 1)
    ϕedge5 = ϕ(x, y, z)[:,:,1]
    # edge 6
    x = range(Ωˣ.min, Ωˣ.max, length = Nx)
    y = range(Ωʸ.min, Ωʸ.max, length = Ny)
    z = range(Ωᶻ.max, Ωᶻ.max, length = 1)
    ϕedge6 = ϕ(x, y, z)[:,:,1]

    file["state"][key] = [ϕedge1, ϕedge2, ϕedge3, ϕedge4, ϕedge5, ϕedge6]
end
