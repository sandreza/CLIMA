#!/usr/bin/env julia --project

include("heat.jl")
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

    ρ = model.ρₒ
    state.ρ = ρ
    state.ρu = ρ * @SVector [-0, -0, -0]
    state.ρθ = ρ * cos(π * z) 

    return nothing
end

#################
# RUN THE TESTS #
#################

vtkpath = abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_box_3D"))

let
    filename  = "heatequation"
    # simulation times
    timeend = FT(0.1) # s
    dt = FT(1e-4) # s
    nout = Int(100)

    # Domain Resolution
    N = 2
    Nˣ = 1
    Nʸ = 1
    Nᶻ = 24

    # Domain size
    Lˣ = 2.0  # m
    Lʸ = 2.0  # m
    Lᶻ = 2.0  # m

    # model params
    cₛ = 0.01 # m/s
    ρₒ = 1 # kg/m³
    μ = 0 # 1e-6,   # m²/s
    ν = 1e-2   # m²/s
    κ = 1e-0   # m²/s
    α = 2e-4   # 1/K
    g = 00.0     # m/s²

    resolution = (; N, Nˣ, Nʸ, Nᶻ)
    domain = (; Lˣ, Lʸ, Lᶻ)
    timespan = (; dt, nout, timeend)
    params = (; cₛ, ρₒ, μ, ν, κ, α, g)

    BC = (
        ClimateMachine.Ocean.OceanBC(
            Impenetrable(NoSlip()),
             TemperatureFlux((state, aux, t) -> ( -κ * 0.0))),
        ClimateMachine.Ocean.OceanBC(
            Impenetrable(KinematicStress(
                (state, aux, t) -> (@SVector [0.00 / state.ρ, -0, -0]),
            )),
            TemperatureFlux((state, aux, t) -> (κ * 0.0)),
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

    f = jldopen(filename * ".jld2", "a+")
    f["grid"] = config.dg.grid
    close(f)

    tic = Base.time()

    run_CNSE(config, resolution, timespan; TimeStepper = SSPRK22Heuns)

    toc = Base.time()
    time = toc - tic
    println(time)
end
