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


include("spatial_model.jl")
include("domain_hooks.jl")
include("CNSE.jl")
include("dissipation_models.jl")

include("abstract_timesteppers.jl")
include("abstract_simulation.jl")

include("./three_dimensional/ThreeDimensionalCompressibleNavierStokesEquations.jl")

ThreeDimensionalCompressibleNavierStokesEquations = ThreeDimensionalCompressibleNavierStokes.CNSE3D

# this does not seem like the right idea
ThreeDimensionalCompressibleNavierStokes.CNSE3D() = ThreeDimensionalCompressibleNavierStokesEquations(1.0,1.0,1.0,1.0,1.0, 1.0, 1.0, 1.0)

# include(pwd() * "/test/Numerics/DGMethods/compressible_navier_stokes_equations/boiler_plate.jl")