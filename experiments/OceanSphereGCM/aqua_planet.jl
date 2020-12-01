using Test
using ClimateMachine
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Grids: polynomialorders
using ClimateMachine.BalanceLaws
using ClimateMachine.Ocean
using ClimateMachine.Ocean.HydrostaticBoussinesq
using ClimateMachine.Ocean.OceanProblems

using CLIMAParameters
using CLIMAParameters.Planet: grav
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

function config_aqua_planet(name, resolution, domain_height, problem, BC)
    problem = problem{FT}(domain_height; BC = BC)

    model = HydrostaticBoussinesqModel{FT}(param_set, problem, cʰ = 1, αᵀ = 0)

    N, Nʰ, Nᶻ = resolution
    resolution = (Nʰ, Nᶻ)

    config = ClimateMachine.OceanSphereGCMConfiguration(
        name,
        N,
        resolution,
        domain_height,
        param_set,
        model,
    )

    return config
end

function run_aqua_planet(driver_config, timespan, Δt; refDat = ())
    grid = driver_config.grid
    vertorder = polynomialorders(grid)[3]
    vert_filter = CutoffFilter(grid, vertorder - 1)
    exp_filter = ExponentialFilter(grid, 1, 8)
    modeldata = (vert_filter = vert_filter, exp_filter = exp_filter)

    solver_type = ClimateMachine.ExplicitSolverType(
        solver_method = LSRK144NiegemannDiehlBusch,
    )

    solver_config = ClimateMachine.SolverConfiguration(
        timespan...,
        driver_config,
        init_on_cpu = true,
        ode_solver_type = solver_type,
        ode_dt = Δt,
        modeldata = modeldata,
    )

    ## Create a callback to report state statistics for main MPIStateArrays
    ## every ntFreq timesteps.
    nt_freq = floor(Int, 1 // 10 * solver_config.timeend / solver_config.dt)
    cb = ClimateMachine.StateCheck.sccreate(
        [
            (solver_config.Q, "3D state"),
            (solver_config.dg.state_auxiliary, "3D aux"),
        ],
        nt_freq;
        prec = 12,
    )

    result = ClimateMachine.invoke!(solver_config; user_callbacks = [cb])

    ## Check results against reference if present
    ClimateMachine.StateCheck.scprintref(cb)
    if length(refDat) > 0
        @test ClimateMachine.StateCheck.scdocheck(cb, refDat)
    end
end
