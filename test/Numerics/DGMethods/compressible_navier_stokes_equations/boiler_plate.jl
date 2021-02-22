using ClimateMachine, MPI
ClimateMachine.init()

using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Grids: polynomialorders
using ClimateMachine.Ocean

using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.BalanceLaws:
    vars_state, Prognostic, Auxiliary, number_states
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.VTK

using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates




include("spatial_model.jl")
include("domain_hooks.jl")
include("CNSE.jl")
include("dissipation_models.jl")

include("abstract_timesteppers.jl")
include("bigfileofstuff.jl")
include("abstract_simulation.jl")

include("./three_dimensional/ThreeDimensionalCompressibleNavierStokesEquations.jl")
# for some reason config needs to be called before anything works
# include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/three_dimensional/run_box.jl")
include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/three_dimensional/box.jl")
Config("", (Nˣ = 1, Nʸ=1, Nᶻ = 1, N = 1), (Lˣ = 1, Lʸ = 1, Lᶻ = 1), (cₛ=1, cᶻ=1, ρₒ=1, μ=0, ν=0, κ=0, α=0, g=0))
# this does not seem like the right idea
ThreeDimensionalCompressibleNavierStokes.CNSE3D() = ThreeDimensionalCompressibleNavierStokes.CNSE3D(1.0,1.0,1.0,1.0,1.0, 1.0, 1.0, 1.0)
ThreeDimensionalCompressibleNavierStokesEquations = ThreeDimensionalCompressibleNavierStokes.CNSE3D

include("simulation_to_clm.jl")
include("simulation_to_run.jl")
