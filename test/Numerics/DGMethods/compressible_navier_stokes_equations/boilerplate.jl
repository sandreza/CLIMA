using MPI
using JLD2
using Test
using Dates
using Printf
using Logging
using StaticArrays
using LinearAlgebra

using ClimateMachine
using ClimateMachine.Mesh.Filters
using ClimateMachine.MPIStateArrays
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Geometry
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws
using ClimateMachine.ODESolvers

# ×(a::SVector, b::SVector) = StaticArrays.cross(a, b)
⋅(a::SVector, b::SVector) = StaticArrays.dot(a, b)
⊗(a::SVector, b::SVector) = a * b'

abstract type AbstractFluid <: BalanceLaw end
struct Fluid <: AbstractFluid end

include("shared_source/domains.jl")
include("shared_source/grids.jl")
include("shared_source/FluidBC.jl")
include("shared_source/abstractions.jl")
include("shared_source/callbacks.jl")

"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
# Description
Gets the (x,y,z) coordinates corresponding to the grid
# Arguments
- `grid`: DiscontinuousSpectralElementGrid
# Return
- `x, y, z`: views of x, y, z coordinates
"""

function evolve!(simulation, spatialmodel; refDat = ())
    Q = simulation.state

    # actually apply initial conditions
    for s in keys(simulation.initial_conditions)
        x, y, z = coordinates(simulation)
        p = spatialmodel.parameters
        ic = simulation.initial_conditions[s]
        ϕ = getproperty(Q, s)
        set_ic!(ϕ, ic, x, y, z, p)
    end

    Ns = polynomialorders(spatialmodel)

    if haskey(spatialmodel.numerics, :overintegration)
        Nover = spatialmodel.numerics.overintegration
    else
        Nover = (0, 0, 0)
    end
    dg = simulation.model

    if sum(Nover) > 0
        cutoff = CutoffFilter(dg.grid, Ns .- (Nover .- 1))
        num_state_prognostic = number_states(dg.balance_law, Prognostic())
        Filters.apply!(Q, 1:num_state_prognostic, dg.grid, cutoff)
    end

    function custom_tendency(tendency, x...; kw...)
        dg(tendency, x...; kw...)
        if sum(Nover) > 0
            cutoff = CutoffFilter(dg.grid, Ns .- (Nover .- 1))
            num_state_prognostic = number_states(dg.balance_law, Prognostic())
            Filters.apply!(tendency, 1:num_state_prognostic, dg.grid, cutoff)
        end
    end

    Δt = simulation.timestepper.timestep
    timestepper = simulation.timestepper.method

    odesolver =
        timestepper(custom_tendency, Q, dt = Δt, t0 = simulation.time.start)

    cbvector = create_callbacks(simulation, odesolver)

    if isempty(cbvector)
        solve!(Q, odesolver; timeend = simulation.time.finish)
    else
        solve!(
            Q,
            odesolver;
            timeend = simulation.time.finish,
            callbacks = cbvector,
        )
    end

    ## Check results against reference
    callbacks = simulation.callbacks
    stcheck = maximum(typeof.(callbacks) .<: StateCheck)
    if stcheck==true
        state_check_ind = argmax(typeof.(callbacks) .<: StateCheck)
        ClimateMachine.StateCheck.scprintref(cbvector[state_check_ind])
        if length(refDat) > 0
            @test ClimateMachine.StateCheck.scdocheck(cbvector[end], refDat)
        end
    end

    return Q
end

function visualize(
    simulation::Simulation;
    statenames = [string(i) for i in 1:size(simulation.state)[2]],
    resolution = (32, 32, 32),
)
    a_, statesize, b_ = size(simulation.state)
    mpistate = simulation.state
    grid = simulation.model.grid
    grid_helper = GridHelper(grid)
    r = coordinates(grid)
    states = []
    ϕ = ScalarField(copy(r[1]), grid_helper)
    r = uniform_grid(Ω, resolution = resolution)
    # statesymbol = vars(Q).names[i] # doesn't work for vectors
    for i in 1:statesize
        ϕ .= mpistate[:, i, :]
        ϕnew = ϕ(r...)
        push!(states, ϕnew)
    end
    visualize([states...], statenames = statenames)
end

# initialized on CPU so not problem, but could do kernel launch?
function set_ic!(ϕ, s::Number, x, y, z, p)
    ϕ .= s
    return nothing
end

function set_ic!(ϕ, s::Function, x, y, z, p)
    a_, states, c_ = size(ϕ)
    @inbounds for i in eachindex(x)
        @inbounds for j in 1:states
            ϕʲ = view(ϕ, :, j, :)
            ϕʲ[i] = s(x[i], y[i], z[i], p)[j]
        end
    end
    return nothing
end
