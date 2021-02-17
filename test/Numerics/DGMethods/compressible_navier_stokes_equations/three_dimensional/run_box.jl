#!/usr/bin/env julia --project

include("box.jl")
ClimateMachine.init()

const FT = Float64

#################
# Initial State #
#################
import ClimateMachine.Ocean: ocean_init_state!

function ocean_init_state!(
    model::ThreeDimensionalCompressibleNavierStokes.CNSE3D,
    state,
    aux,
    localgeo,
    t,
)

    x = aux.x
    y = aux.y
    z = aux.z
    cₛ = model.cₛ
    ρ = model.ρₒ * ( 1 +  ( 2e-3 / cₛ^2) * z^2 / 200.0 )
    # 0.01 * z + 0.01 * 2e-3 / cₛ^2 z^3 / 200
    # 0.01 * (1 + 0.01 * 2e-3 / cₛ^2 * 3 * 100^2 / 200)
    state.ρ = ρ
    ϵ = 1e-6
    state.ρu = ρ * @SVector [ϵ*sin(2π*y/100)*sin(2π*z/100), ϵ*sin(2π*x/100)*sin(2π*z/100), -0]
    state.ρθ = ρ * (0.01 * z )

    return nothing
end

#################
# RUN THE TESTS #
#################

vtkpath = abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_box_3D"))

let
    cₛ = 0.2 # sqrt(1) # m/s # dont forget to change initial condition
    filename  = "cool_the_box_3"
    # simulation times
    timeend = FT(4 * 24 * 60 * 60) # s
    dt = FT(0.125 / cₛ ) # s # 0.125 ≈ Lˣ / ( Nˣ * (N + 1)^2) * 0.3
    nout = round(Int, 60 * 60 * 4 / dt) # output every 4 hours or so

    # Domain Resolution
    N = 3
    Nˣ = 8*2
    Nʸ = 8*2
    Nᶻ = 8*2

    # Domain size
    Lˣ = 100.0  # m
    Lʸ = 100.0  # m
    Lᶻ = 100.0  # m

    # model params
    ρₒ = 1 # kg/m³
    μ = 0 # 1e-6,   # m²/s
    ν = 5e-3   # m²/s ν = 1e-2, κ = 1e-4 worked
    κ = 1e-4   # m²/s
    α = 2e-4   # 1/K
    g = 10.0   # m/s²

    resolution = (; N, Nˣ, Nʸ, Nᶻ)
    domain = (; Lˣ, Lʸ, Lᶻ)
    timespan = (; dt, nout, timeend)
    params = (; cₛ, ρₒ, μ, ν, κ, α, g)

    BC = (
        ClimateMachine.Ocean.OceanBC(
            Impenetrable(NoSlip()),
             TemperatureFlux((state, aux, t) -> ( (0.01 + 0.01 * 2e-3 / cₛ^2 * 3 * 100^2 / 200) * κ))),
        ClimateMachine.Ocean.OceanBC(
            Impenetrable(KinematicStress(
                (state, aux, t) -> (@SVector [0.00 / state.ρ, -0, -0]),
            )),
            TemperatureFlux((state, aux, t) -> (1e-5)),
        ),
    )

    config = Config(
        filename,
        resolution,
        domain,
        params;
        numerical_flux_first_order = RoeNumericalFlux(),
        Nover = 1,
        periodicity = (true, true, false),
        boundary = ((0, 0), (0, 0), (1, 2)),
        boundary_conditons = BC,
    )

    # f = jldopen(filename * ".jld2", "a+")
    # f["grid"] = config.dg.grid
    # close(f)

    tic = Base.time()

    run_CNSE(config, resolution, timespan; TimeStepper = SSPRK22Heuns)

    toc = Base.time()
    time = toc - tic
    println(time)
end
