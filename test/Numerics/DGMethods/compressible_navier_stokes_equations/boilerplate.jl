using MPI
using JLD2
using Test
using Dates
using Printf
using Logging
using StaticArrays
using LinearAlgebra

using ClimateMachine
using ClimateMachine.MPIStateArrays
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Geometry
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws
using ClimateMachine.ODESolvers

import ClimateMachine.Mesh.Grids: polynomialorders

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
function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end

function evolve!(simulation, spatialmodel)
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

    odesolver = timestepper(
        custom_tendency,
        Q,
        dt = Δt,
        t0 = simulation.simulation_time[1],
    )

    cbvector = [nothing] # create_callbacks(simulation)

    if cbvector == [nothing]
        solve!(Q, odesolver; timeend = simulation.simulation_time[2])
    else
        solve!(
            Q,
            odesolver;
            timeend = simulation.simulation_time[2],
            callbacks = cbvector,
        )
    end
    return Q
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
