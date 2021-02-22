using Test
using StaticArrays
using LinearAlgebra

using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Geometry
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws
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
import ClimateMachine.NumericalFluxes: numerical_flux_first_order!

×(a::SVector, b::SVector) = StaticArrays.cross(a, b)
⋅(a::SVector, b::SVector) = StaticArrays.dot(a, b)
⊗(a::SVector, b::SVector) = a * b'