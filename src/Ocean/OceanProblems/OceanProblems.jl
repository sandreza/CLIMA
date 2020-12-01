module OceanProblems

export SimpleBox, Fixed, Rotating, HomogeneousBox, OceanGyre
export SimpleSphere

using StaticArrays
using CLIMAParameters.Planet: grav

using ...Problems

using ..Ocean
using ..HydrostaticBoussinesq
using ..ShallowWater
using ..SplitExplicit01

using CLIMAParameters
using CLIMAParameters.Planet: planet_radius

import ..Ocean:
    ocean_init_state!,
    ocean_init_aux!,
    kinematic_stress,
    surface_flux,
    coriolis_parameter

HBModel = HydrostaticBoussinesqModel
SWModel = ShallowWaterModel

"""
    ocean_init_aux!(::HBModel, ::AbstractOceanProblem)

save y coordinate for computing coriolis, wind stress, and sea surface temperature

# Arguments
- `m`: model object to dispatch on and get viscosities and diffusivities
- `p`: problem object to dispatch on and get additional parameters
- `A`: auxiliary state vector
- `geom`: geometry stuff
"""
function ocean_init_aux!(::HBModel, ::AbstractOceanProblem, A, geom)
    FT = eltype(A)
    @inbounds A.y = geom.coord[2]

    # needed for proper CFL condition calculation
    A.w = -0
    A.pkin = -0
    A.wz0 = -0

    A.uᵈ = @SVector [-0, -0]
    A.ΔGᵘ = @SVector [-0, -0]

    return nothing
end

function ocean_init_aux!(::OceanModel, ::AbstractOceanProblem, A, geom)
    FT = eltype(A)
    @inbounds A.y = geom.coord[2]

    # needed for proper CFL condition calculation
    A.w = -0
    A.pkin = -0
    A.wz0 = -0

    A.u_d = @SVector [-0, -0]
    A.ΔGu = @SVector [-0, -0]

    return nothing
end

function ocean_init_aux!(::SWModel, ::AbstractOceanProblem, A, geom)
    @inbounds A.y = geom.coord[2]

    A.Gᵁ = @SVector [-0, -0]
    A.Δu = @SVector [-0, -0]

    return nothing
end

function ocean_init_aux!(::BarotropicModel, ::AbstractOceanProblem, A, geom)
    @inbounds A.y = geom.coord[2]

    A.Gᵁ = @SVector [-0, -0]
    A.U_c = @SVector [-0, -0]
    A.η_c = -0
    A.U_s = @SVector [-0, -0]
    A.η_s = -0
    A.Δu = @SVector [-0, -0]
    A.η_diag = -0
    A.Δη = -0

    return nothing
end

@inline coriolis_parameter(
    m::Union{HBModel, SWModel, OceanModel, BarotropicModel},
    ::AbstractOceanProblem,
    y,
) = -0

include("SimpleBoxProblem.jl")
include("SimpleSphereProblem.jl")

end
