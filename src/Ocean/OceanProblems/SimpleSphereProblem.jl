abstract type AbstractSimpleSphereProblem <: AbstractOceanProblem end

struct SimpleSphere{T, BC} <: AbstractSimpleSphereProblem
    H::T
    boundary_conditions::BC
    function SimpleSphere{FT}(
        H; # m
        BC = (
            OceanBC(Impenetrable(FreeSlip()), Insulating()),
            OceanBC(Penetrable(FreeSlip()), Insulating()),
        ),
    ) where {FT <: AbstractFloat}
        return new{FT, typeof(BC)}(H, BC)
    end
end

function ocean_init_state!(m::HBModel, p::SimpleSphere, Q, A, localgeo, t)
    Q.u = @SVector [-0, -0]
    Q.η = -0

    @inbounds z = localgeo.coord[3] - planet_radius(m.param_set)

    Q.θ = 20 # * (1 + z / p.H)

    return nothing
end

function ocean_init_state!(::SWModel, ::SimpleSphere, Q, A, localgeo, t)
    Q.U = @SVector [-0, -0]
    Q.η = -0

    return nothing
end
