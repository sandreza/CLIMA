# need to make the following
# to check for tmorrow
# make sure that specifying boundaries strangely doesn't change the solution
#=

ρu_bcs = (bottom = NoSlip(), top = FreeSlip())
# temperature, default no-flux
ρθ_bcs = (bottom = Insulating(), top = TemperatureFlux())
# combine
boundary_conditions = (ρθ = ρθ_bcs, ρu = ρu_bcs)

BC = (
    ClimateMachine.Ocean.OceanBC(Impenetrable(NoSlip()), Insulating()),
    ClimateMachine.Ocean.OceanBC(
        Impenetrable(KinematicStress(
            (state, aux, t) -> (@SVector [0.01 / state.ρ, -0, -0]),
        )),
        TemperatureFlux((state, aux, t) -> (0.1)),
    ),
)
((0, 0), (0, 0), (1, 2)),
# Redundancies, Flux
=#

using Test
using StaticArrays
using LinearAlgebra

using ClimateMachine.Ocean
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Geometry
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws
using ClimateMachine.Orientations
using ClimateMachine.Ocean: coriolis_parameter
using ClimateMachine.Mesh.Geometry: LocalGeometry
using ClimateMachine.MPIStateArrays: MPIStateArray

import ClimateMachine.BalanceLaws:
    vars_state,
    init_state_prognostic!,
    init_state_auxiliary!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    flux_first_order!,
    flux_second_order!,
    source!,
    wavespeed,
    boundary_conditions,
    boundary_state!
import ClimateMachine.Ocean:
    ocean_init_state!,
    ocean_init_aux!,
    ocean_boundary_state!,
    _ocean_boundary_state!
import ClimateMachine.NumericalFluxes: numerical_flux_first_order!
import ClimateMachine.Orientations: vertical_unit_vector





#=
function convention(a::NamedTuple{(:vertical, :horizontal), T}, ::Val{3}) where T
    return (a.horizontal, a.horizontal, a.vertical)
end

function convention(a::Number, ::Val{3})
    return (a, a, a)
end

function convention(a::NamedTuple{(:vertical, :horizontal), T}, ::Val{2}) where T
    return (a.horizontal, a.vertical)
end

function convention(a::Number, ::Val{2})
    return (a, a)
end

function convention(a::Tuple)
    return a
end




elements = (vertical = 4, horizontal = 8)
convention(elements, Val(3))
convention(elements, Val(2))
=#
