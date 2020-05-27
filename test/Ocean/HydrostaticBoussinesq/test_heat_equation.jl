using Test
using CLIMA
using CLIMA.HydrostaticBoussinesq
using CLIMA.GenericCallbacks
using CLIMA.ODESolvers
using CLIMA.Mesh.Filters
using CLIMA.PlanetParameters
using CLIMA.VariableTemplates
using CLIMA.Mesh.Grids: polynomialorder

function config_simple_box(FT, N, resolution, dimensions)
    prob = HeatedBox{FT}(dimensions...)

    κ = 0.001
    model =
        HydrostaticBoussinesqModel{FT}(prob, αᵀ = 0, νʰ = 5e-3, κʰ = κ, κᶻ = κ)

    config =
        CLIMA.OceanBoxGCMConfiguration("heat_equation", N, resolution, model)

    return config
end

function main(; imex::Bool = false, Δt = 60)
    CLIMA.init()

    FT = Float64

    # DG polynomial order
    N = Int(4)

    # Domain resolution and size
    Nˣ = Int(2)
    Nʸ = Int(2)
    Nᶻ = Int(12)
    resolution = (Nˣ, Nʸ, Nᶻ)

    Lˣ = 10    # m
    Lʸ = 10    # m
    H = 60    # m
    dimensions = (Lˣ, Lʸ, H)

    timestart = FT(0)    # s
    timeout = FT(60)   # s
    timeend = FT(3600) # s

    if imex
        solver_type = CLIMA.IMEXSolverType(linear_model = LinearHBModel)
    else
        solver_type =
            CLIMA.ExplicitSolverType(solver_method = LSRK144NiegemannDiehlBusch)
    end

    driver_config = config_simple_box(FT, N, resolution, dimensions)

    grid = driver_config.grid
    vert_filter = CutoffFilter(grid, polynomialorder(grid) - 1)
    exp_filter = ExponentialFilter(grid, 1, 8)
    modeldata = (vert_filter = vert_filter, exp_filter = exp_filter)

    solver_config = CLIMA.setup_solver(
        timestart,
        timeend,
        driver_config,
        init_on_cpu = true,
        ode_solver_type = solver_type,
        # ode_dt = FT(Δt),
        modeldata = modeldata,
    )

    CLIMA.Settings.enable_vtk = true
    CLIMA.Settings.vtk_interval = ceil(timeout / solver_config.dt)

    result = CLIMA.invoke!(solver_config)

    @test true
end

@testset "$(@__FILE__)" begin
    main(imex = false, Δt = 6)
    # main(imex=true)
end
