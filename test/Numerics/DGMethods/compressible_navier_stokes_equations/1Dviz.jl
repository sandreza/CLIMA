using ClimateMachine, JLD2, GLMakie

filename = "heatequation"
filename = "cool_the_box"
f = jldopen(filename * ".jld2", "r+")

include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/bigfileofstuff.jl")
include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/ScalarFields.jl")

nout = 60 + 1
dg_grid = f["grid"]
gridhelper = GridHelper(dg_grid)   
x, y, z = coordinates(dg_grid)
ϕ =  ScalarField(copy(x), gridhelper)
#=
# analytic solution for heat equation
dt = 1e-4 * 100
analytic = zeros(size(x)..., nout)
for i in 1:nout
    # @. analytic[:,:,i] = 
    tmp = @. cos(π*z) * exp(-(π)^2 * (i-1) * dt)
    Q = f[string(i-1)][:,5,:]
    println(norm(tmp-Q)/norm(Q))
end
i = 10
Q = f[string(i-1)][:,5,:]
tmp = @. cos(π*z) * exp(-(π)^2 * (i-1) * dt)
norm(Q - tmp) / norm(Q)
=#
##
Ω = (extrema(x), extrema(y), extrema(z))
newx = range(Ω[1][1], Ω[1][2], length = 2 )
newy = range(Ω[2][1], Ω[2][2], length = 2 )
newz = range(Ω[3][1], Ω[3][2], length = 32 )
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

##
# oldρθ = ρθ
fig = Figure()
ax1 = fig[2:4, 2:4] = Axis(fig, title = "Heat Equation", xlabel = "Temperature", ylabel = "Depth")

timeslider = Slider(fig, range = Int.(collect(1:nout)))
timenode = timeslider.value

state = @lift(ρθ[1,1,:,$timenode])
tmp = lines!(ax1, state, -100..0, color = :red, linewidth = 5)

#=
oldstate = @lift(oldρθ[1,1,:,$timenode])
tmp2 = lines!(ax1, oldstate, -100..0, color = :blue, linewidth = 5)

fig[2, 1] = vgrid!(
        Label(fig, "Time", width = nothing),
        timeslider,
        Legend(fig,
        [tmp, tmp2],
        ["Overintegration", "Underintegration",])
)
=#
fig[2, 1] = vgrid!(
        Label(fig, "Time", width = nothing),
        timeslider,
)
display(fig)
##
record(fig, "heat.mp4", 1:nout, framerate=10) do n
    timenode[] = n
end

##
# Target Design Template for Idealized Test Cases
# domain / grid
Ω = Atmosphere(height = 30.0km, radius = 6378km)
grid = DiscontinuousSpectralElementGrid(
    domain = Ω,
    elements = (vertical = 4, horizontal = 8)
    polynomialorder = (vertical = 1, horizontal = 5)
)

##


