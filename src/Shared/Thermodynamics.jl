module Thermodynamics
"""
    Thermodynamics

Module containing various thermodynamic functions, e.g., for air pressure
(atmosphere equation of state), latent heats of phase transitions, and
saturation vapor pressures.
"""

# Tests to add:
# air_pressure([1, 1], [1, 1], [1, 0], [0, 0]) = [R_v, R_d]
# latent_heat_vapor(T_0) = L_v0
# latent_heat_fusion(T_0) = L_f0
# latent_heat_sublim(T_0) = L_s0
# sat_vapor_press_liquid(T_triple) = sat_vapor_press_triple
# sat_vapor_press_ice(T_triple) = sat_vapor_press_triple

using PlanetParameters

# Atmospheric equation of state
export air_pressure

# Latent heats
export latent_heat_vapor, latent_heat_sublim, latent_heat_fusion

# Saturation vapor pressures over liquid and ice
export sat_vapor_press_liquid, sat_vapor_press_ice

function air_pressure(T, density, q_t, q_c)
    """
        air_pressure(T, density, q_t, q_c)

    Computes air pressure from the equation of state (ideal gas law) given
    the temperature T, density, total water specific humidity q_t, and condensed
    water vapor specific humidity q_c (so that the water vapor specific humidity
    is q_v = q_t - q_c).

    """

    p = R_d * density .* T .*( 1.0 .+  (molmass_ratio - 1.0) * q_t .- molmass_ratio * q_c )

end

function latent_heat_vapor(T)

    L_v = latent_heat_generic(T, L_v0, cp_v - cp_l)

end

function latent_heat_sublim(T)

    L_s = latent_heat_generic(T, L_s0, cp_v - cp_i)

end

function latent_heat_fusion(T)

    L_f = latent_heat_generic(T, L_f0, cp_l - cp_i)

end

function latent_heat_generic(T, L_0, cp_diff)
    """
        latent_heat_generic(T, L_0, cp_diff)

    Computes the latent heat of a generic phase transition between two
    phases using Kirchhoff's relation. L_0 is the latent heat of the phase
    transition at T_0, and cp_diff is the difference between the isobaric
    specific heat capacities (cp in higher-temperature phase minus cp in
    lower-temperature phase). The isobaric specific heat capacities are
    assumed to be constant.
    """

    L = L_0 .+ cp_diff * (T .- T_0)

end

function sat_vapor_press_liquid(T)

    es = sat_vapor_press_generic(T, L_v0, cp_v - cp_l)

end

function sat_vapor_press_ice(T)

    es = sat_vapor_press_generic(T, L_s0, cp_v - cp_i)

end

function sat_vapor_press_generic(T, L_0, cp_diff)
    """
        sat_vapor_press_generic(T, L_0, cp_diff)

    Computes vapor pressure by integration of Clausius-Clepeyron equation
    from the triple point T_triple, using Kirchhoff's relation

        L = L_0 + cp_diff * (T - T_0)

    for the latent heat L and assuming constant isobaric specific heats of
    the phases. This gives a linear dependence of latent heat on temperature
    around a reference temperature T_0, which can be integrated analytically.
    That is, the saturation vapor pressure is the solution of

        dlog(e_s)/dT = [L_0 + cp_diff(T-T_0)]/(R_v*T^2)

    """

    es = sat_vapor_press_triple * exp.(
        (L_0 - cp_diff*T_0)/R_v * (1.0/T_triple .- 1.0./T)
        .+ cp_diff/R_v*log.(T/T_triple)
        )

end

end