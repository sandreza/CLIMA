

##
# surface 
Ωˣ = cpu_grid.domain[1]
Ωʸ = cpu_grid.domain[2]
Ωᶻ = cpu_grid.domain[3]
xsurf = range(Ωˣ.min, Ωˣ.max, length = 64)
ysurf = range(Ωʸ.min, Ωʸ.max, length = 64)
zsurf = range(Ωᶻ.max, Ωᶻ.max, length = 1)
ϕsurf = ϕ(xsurf, ysurf, zsurf)
clims = extrema(ϕsurf)
zscale = 100
fig = Figure(resolution = (1920, 1080))
ax = fig[1,1] = LScene(fig, title= "Baroclinic Adjustment")


# edge 1
x = range(Ωˣ.min, Ωˣ.max, length = 64)
y = range(Ωʸ.min, Ωʸ.min, length = 1)
z = range(Ωᶻ.min, Ωᶻ.max, length = 16)
ϕedge1 = ϕ(x, y, z)[:,1,:]
surface!(ax, x, z .* zscale, ϕedge1, transformation = (:xz, Ωʸ.min),  colorrange = clims, colormap = :balance, show_axis=false)

# edge 2
x = range(Ωˣ.min, Ωˣ.max, length = 64)
y = range(Ωʸ.min, Ωʸ.max, length = 64)
z = range(Ωᶻ.min, Ωᶻ.max, length = 16)
ϕedge2 = ϕ(x, y, z)[:,1,:]
surface!(ax, x, z .* zscale, ϕedge2, transformation = (:xz, Ωʸ.max),  colorrange = clims, colormap = :balance)


# edge 3
x = range(Ωˣ.min, Ωˣ.min, length = 1)
y = range(Ωʸ.min, Ωʸ.max, length = 64)
z = range(Ωᶻ.min, Ωᶻ.max, length = 16)
ϕedge3 = ϕ(x, y, z)[1,:,:]
surface!(ax, y, z .* zscale, ϕedge3, transformation = (:yz, Ωˣ.min),  colorrange = clims, colormap = :balance)

# edge 4
x = range(Ωˣ.max, Ωˣ.max, length = 1)
y = range(Ωʸ.min, Ωʸ.max, length = 64)
z = range(Ωᶻ.min, Ωᶻ.max, length = 16)
ϕedge4 = ϕ(x, y, z)[1,:,:]
surface!(ax, y, z .* zscale, ϕedge4, transformation = (:yz, Ωˣ.max),  colorrange = clims, colormap = :balance)

# edge 5
x = range(Ωˣ.min, Ωˣ.max, length = 64)
y = range(Ωʸ.min, Ωʸ.max, length = 64)
z = range(Ωᶻ.min, Ωᶻ.min, length = 1)
ϕedge5 = ϕ(x, y, z)[:,:,1]
surface!(ax, x, y, ϕedge5, transformation = (:xy, Ωᶻ.min *  zscale), colorrange = clims, colormap = :balance)


# edge 6
x = range(Ωˣ.min, Ωˣ.max, length = 64)
y = range(Ωʸ.min, Ωʸ.max, length = 64)
z = range(Ωᶻ.max, Ωᶻ.max, length = 1)
ϕedge6 = ϕ(x, y, z)[:,:,1]
surface!(ax, x, y, ϕedge6, transformation = (:xy, Ωᶻ.max *  zscale), colorrange = clims, colormap = :balance)

##
seconds = 5
fps = 30
frames = round(Int, fps * seconds )
GLMakie.record(fig.scene, pwd() * "/example.mp4"; framerate = fps) do io
    for i = 1:frames
        sleep(1/fps)
        recordframe!(io)
    end
end
