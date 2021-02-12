using ClimateMachine, JLD2, GLMakie

filename = "cool_the_box"
f = jldopen(filename * ".jld2", "r+")

include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/bigfileofstuff.jl")
include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/ScalarFields.jl")

nout = 20 + 1
dg_grid = f["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
xC, yC, zC = cellcenters(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)

Ω = (extrema(x), extrema(y), extrema(z))
newx = range(Ω[1][1], Ω[1][2], length = 6 )
newy = range(Ω[2][1], Ω[2][2], length = 6 )
newz = range(Ω[3][1], Ω[3][2], length = 6 )
##
ρ  = zeros(length(newx), length(newy), length(newz), nout)
ρu = zeros(length(newx), length(newy), length(newz), nout)
ρv = zeros(length(newx), length(newy), length(newz), nout)
ρw = zeros(length(newx), length(newy), length(newz), nout)
ρθ = zeros(length(newx), length(newy), length(newz), nout)
tic = time()
for i in 0:(nout-1)
    Q = f[string(i)]
    ϕ .= Q[:,1,:]
    ρ[:,:,:, i+1]  = ϕ(newx, newy, newz, threads = true)
    ϕ .= Q[:,2,:]
    ρu[:,:,:, i+1] = ϕ(newx, newy, newz, threads = true)
    ϕ .= Q[:,3,:]
    ρv[:,:,:, i+1] = ϕ(newx, newy, newz, threads = true)
    ϕ .= Q[:,4,:]
    ρw[:,:,:, i+1] = ϕ(newx, newy, newz, threads = true)
    ϕ .= Q[:,5,:]
    ρθ[:,:,:, i+1] = ϕ(newx, newy, newz, threads = true)
end
toc = time()
close(f)
println("time to interpolate is $(toc-tic)")

limits = FRect(19, -100, 1, 100)
scene = lines(ρθ[1,1,:,1], -100..0, limits = limits)
axis = scene.axis # get the axis object from the scene
axis.xlabel = "Temperature"
axis.ylabel = "Depth"
display(scene)