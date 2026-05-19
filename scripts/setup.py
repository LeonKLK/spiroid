"""
Initial conditions for input into spiroid.

A simulation input config file is generated for all combinations
of the values of the star, planet and effects.

"""

##############################################################
########################!!! WARNING !!!#######################
##############################################################
"""
Do _NOT_ change the order of dictionary keys in this file.
The values are unpacked in an order dependent way in config.py
to create the initial conditions.
"""
##############################################################
########################!!! WARNING !!!#######################
##############################################################


import math
import sys

sys.dont_write_bytecode = True
from spiroid.configmaker import make_configs
from units import (
    AU,
    SECONDS_IN_YEAR,
    SOLAR_MASS,
    SOLAR_RADIUS,
    JUPITER_MASS,
    JUPITER_RADIUS,
    EARTH_MASS,
    EARTH_RADIUS,
)


def simulator_setup():
    ##############################################################
    ####################### SIMULATOR SETUP ######################
    ##############################################################

    simulation = {
        # The prefix simulation name.
        "name": "1Msun_1.91d_a03",
        # Decription of the science case.
        "decription": "",
        # Simulation start time, seconds (from years).
        "start_time": SECONDS_IN_YEAR * 5e6,
        # Simulation end time, seconds (from years).
        "final_time": SECONDS_IN_YEAR * 5e7,
    }

    # seconds (from years) — set before start so disk is already dissipated at t=5 Myr
    disk_lifetime = SECONDS_IN_YEAR * 1e6

    return (simulation, disk_lifetime)


def effect_setup():
    # Enables or disables certain effects for all simulations.
    # Must be [True], [False] or [True, False].
    effects = {
        "MAGNETIC_EFFECT_ENABLED": [False],
        "STAR_EVOLUTION_ENABLED": [True],
        # Kaula stellar tide
        "STAR_TIDES_ENABLED": [True],
        "STAR_TIDES_MODEL": "KaulaTides",
        # Kaula planetary tides
        "PLANET_TIDES_ENABLED": [False],
        # Disable wind for testing conservation of angular momentum
        "WIND_ENABLED": [True],
    }

    return effects


def planet_setup(effects):
    ##############################################################
    ####################### PLANET SETUP #########################
    ##############################################################
    planet_base = {
        # kg
        "mass": [5 * (1./332946.050895) * 1.9818e30],
        # m
        "radius": [1.0014e+7],
        # m (from AU)
        "semi_major_axis": [AU * x for x in [0.03]],
        "magnetic_field": [None],  # Do not edit.
    }

    if effects["MAGNETIC_EFFECT_ENABLED"]:
        # Gauss
        planet_base["magnetic_field"] = [10.0]

    if effects["PLANET_TIDES_ENABLED"]:
        # For Kaula
        planet_base.update(
            {
                # rad.s
                "spin": [8.093879511357418e-07],
                # No units
                "eccentricity": [0.005],
                # rad
                "inclination": [0.3490658503988659],
                # rad
                "longitude_ascending_node": [1.0],
                # rad
                "pericentre_omega": [0.0],
                # rad
                "spin_inclination": [0.34906584951436426],
                # No units
                "radius_of_gyration": [0.33070368308499226],
                "type_and_file": [
                    (
                        "Solid",
                        "examples/data/planet/tides/kaula/leconte2015_steinberger.csv",
                    )
                ],
            }
        )
    elif effects["STAR_TIDES_ENABLED"] and effects.get("STAR_TIDES_MODEL") == "KaulaTides":
        # Star Kaula tide requires orbit geometry even though the planet itself is not deformed.
        planet_base.update(
            {
                "eccentricity": [1e-6],
                "inclination": [0.0],
                "longitude_ascending_node": [0.0],
                "pericentre_omega": [0.0],
                "spin_inclination": [0.0],
            }
        )

    return planet_base


def star_setup(effects):
    ##############################################################
    ####################### STAR SETUP ###########################
    ##############################################################
    star_base = {
        "mass": [None],  # Do not edit.
        "radius": [None],  # Do not edit.
        # rad.s-1  (1.91 day rotation period)
        "spin": [2 * math.pi / (1.91 * 86400)],
        # seconds (from years)
        "core_envelope_coupling_constant": [SECONDS_IN_YEAR * x for x in [1.171e7]],
        "footpoint_conductance": [None],  # Do not edit.
        "evolution": [None],  # Do not edit.
        "sigma_bar": [None],  # Do not edit.
    }

    if effects["MAGNETIC_EFFECT_ENABLED"]:
        # Ohm-1
        star_base["footpoint_conductance"] = [5.8e4]

    if effects["STAR_EVOLUTION_ENABLED"]:
        star_base["evolution"] = [
            {"Starevol": {"star_file_path": "examples/data/star/evolution/savgol_10.csv"}},
        ]
    else:
        # Set the initial star values that would otherwise be provided by savgol/mesa data if evolution were enabled.
        # Must be non-zero (to avoid NaN).

        # kg (from Msun)
        star_base["mass"] = [SOLAR_MASS * x for x in [0.8]]
        # m (from SOLAR_RADIUS)
        star_base["radius"] = [SOLAR_RADIUS * x for x in [1.0]]
        # No units
        star_base["radiative_moment_of_inertia"] = [1.0]
        star_base["convective_moment_of_inertia"] = [1.0]

    if effects["STAR_TIDES_ENABLED"] and effects.get("STAR_TIDES_MODEL") == "KaulaTides":
        star_base["stellar_convective_file"] = [
            "examples/data/star/tides/kaula/alpha0.51600_P1p2_Ek1p5em6.json"
        ]
        # rad.s-1 — reference spin rate at which the love number spectrum was computed (1.2 day period)
        star_base["stellar_spectrum_spin_rate"] = [2 * math.pi / (1.2 * 86400)]
    elif effects["STAR_TIDES_ENABLED"]:
        star_base["sigma_bar"] = [1.0e-6]

    return star_base


def integrator_setup():
    ##############################################################
    #################### INTEGRATOR SETUP ########################
    ##############################################################

    # Output filtered each timestep that occurs greater than `absolute_tolerance`
    # or if the percentage difference between two timesteps exceeds `relative_tolerance`,
    # in which case the `relative_tolerance` is increased by the additive `scaling_factor`
    # Used to simulate logscale output.
    filter = {
        "Filtered": {
            "absolute_tolerance": SECONDS_IN_YEAR * 1e3,
            "relative_tolerance": 1.15,
            "incremental_scaling_factor": 0.002,
            "decremental_scaling_factor": 0.0,
        }
    }

    odex = {
        "Odex": {
            "step_size_reduction_factor": 0.66666666666666666666,
            "step_size_selection_b": 2.0,
            "step_size_max": SECONDS_IN_YEAR * 5e5,
            "max_integration_steps": 100000000,
            # Uncomment the entire solution_output for Dense output
            "solution_output": filter,
        }
    }

    # Kaula tides on the planet modifies the integrator defaults to avoid timestep issue.
    odex_kaula = {
        "Odex": {
            "step_size_reduction_factor": 0.66666666666666666666,
            "step_size_selection_b": 2.0,
            "step_size_max": SECONDS_IN_YEAR * 1e3,
            "max_integration_steps": 500000000,
            "step_control_safety_a": 0.05,
            "step_control_safety_b": 0.2,
        }
    }

    dopri853 = {
        "Dopri853": {
            "step_controller": {
                "relative_tolerance": 1e-16,
                "absolute_tolerance": 1e-16,
                "step_size_factor_min": 0.3333333333333333,
                "step_size_factor_max": 6.0,
                "step_size_error_factor": 0.9,
                "step_size_max": SECONDS_IN_YEAR * 1e3,
                "alpha": 0.125,
                "beta": 0.0,
            },
            "step_size_underflow": None,
            "stiffness_test": "Disabled",
            "max_integration_steps": 500000000,
            # Uncomment the entire solution_output for Dense output
            "solution_output": filter,
        }
    }

    # Uncomment only the desired integrator.
    # return odex
    # return odex_kaula
    return dopri853


if __name__ == "__main__":
    make_configs(simulator_setup, effect_setup, planet_setup, star_setup, integrator_setup)
