"""
Initial conditions for input into spiroid: stellar dynamical tide (kaula model) on a rotator.

A 1 Msun star started as a fast, median or slow rotator (Gallet & Bouvier 2013 or 2015 initial
conditions, see scripts/setup_fast_rotator.py and rotators.md), with stellar evolution, wind and
the stellar dynamical tide read from one of the tabulated love number spectra in
examples/data/star/tides/kaula/ (named by age, alpha, beta and rotation period). The orbiting
body is a Jupiter-mass planet on a circular orbit, with no tide of its own.

Pick the rotator (ROTATOR, REFERENCE) and the spectrum (SPECTRUM_FILE) below.

usage: python3 stellar_kaula_setup.py path/to/output/folder        (run from scripts/)
       spiroid -b --output-format jsonl path/to/output/folder path/to/results   (from the repo root)

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
)


TWO_PI = 6.283185307179586

# Rotator classes for a 1 Msun star: initial period (d), disk-locking lifetime (yr) and
# core-envelope coupling time (yr).
#   GB13: Gallet & Bouvier 2013, A&A 556, A36.
#   GB15: Gallet & Bouvier 2015, A&A 577, A98, Table 2.
ROTATORS = {
    "GB13": {
        "slow": {"period_d": 10.0, "disk_lifetime_yr": 5.0e6, "coupling_yr": 30e6},
        "median": {"period_d": 7.0, "disk_lifetime_yr": 5.0e6, "coupling_yr": 28e6},
        "fast": {"period_d": 1.4, "disk_lifetime_yr": 2.5e6, "coupling_yr": 12e6},
    },
    "GB15": {
        "slow": {"period_d": 8.0, "disk_lifetime_yr": 9.0e6, "coupling_yr": 30e6},
        "median": {"period_d": 5.0, "disk_lifetime_yr": 6.0e6, "coupling_yr": 28e6},
        "fast": {"period_d": 1.4, "disk_lifetime_yr": 2.0e6, "coupling_yr": 10e6},
    },
}
REFERENCE = "GB15"  # "GB13" or "GB15"
ROTATOR = "fast"  # "slow", "median" or "fast"

# Stellar track, turnover-time option and the wind prefactor (J) calibrated with it
# (solar period 25.38 d at 4.567 Gyr).
TRACK = {
    "prefactor": 8e23,
    "evolution": {"Starevol": {"star_file_path": "examples/data/star/evolution/starevol_m1p00_2026.csv",
                               "convective_turnover_time_mode": "FromFile"}},
}

# Love number spectrum of the stellar dynamical tide (inertial waves in the convective
# envelope), see examples/data/star/tides/kaula/. The file carries the stellar spin it was
# computed at; spiroid rescales it to the current spin during the run.
SPECTRUM_FILE = "examples/data/star/tides/kaula/spectrum_age41.5Myr_alpha0.74610_beta0.98032_P0.83820006d.json"


def simulator_setup():
    ##############################################################
    ####################### SIMULATOR SETUP ######################
    ##############################################################

    simulation = {
        # The prefix simulation name.
        "name": f"stellar_kaula_{ROTATOR}_{REFERENCE.lower()}",
        # Decription of the science case.
        "decription": "",
        # Simulation start time, seconds (from years).
        "start_time": SECONDS_IN_YEAR * 2.0e6,  # the FromFile mode needs the radiative core (1.69 Myr on this track); earlier starts stall
        # Simulation end time, seconds (from years).
        "final_time": SECONDS_IN_YEAR * 5.0e9,
    }

    # seconds (from years)
    disk_lifetime = SECONDS_IN_YEAR * ROTATORS[REFERENCE][ROTATOR]["disk_lifetime_yr"]

    return (simulation, disk_lifetime)


def effect_setup():
    # Enables or disables certain effects for all simulations.
    # Must be [True], [False] or [True, False].
    effects = {
        "MAGNETIC_EFFECT_ENABLED": [False],
        "STAR_EVOLUTION_ENABLED": [True],
        # Stellar tides, one model at most: the kaula stellar tide reads SPECTRUM_FILE above.
        "STAR_TIDES_KAULA_ENABLED": [True],
        "STAR_TIDES_CTL_ENABLED": [False],
        # Planetary tides, one model at most.
        "PLANET_TIDES_KAULA_ENABLED": [False],
        "PLANET_TIDES_CTL_ENABLED": [False],
        # Disable wind for testing conservation of angular momentum
        "WIND_ENABLED": [True],
        "GENERAL_RELATIVITY_ENABLED": [False],
    }

    return effects


def planet_setup(effects):
    ##############################################################
    ####################### PLANET SETUP #########################
    ##############################################################
    planet_base = {
        # kg
        "mass": [JUPITER_MASS],
        # m
        "radius": [JUPITER_RADIUS],
        # m (from AU)
        "semi_major_axis": [AU * x for x in [0.03]],
        "magnetic_field": [None],  # Do not edit.
    }

    if effects["MAGNETIC_EFFECT_ENABLED"]:
        # Gauss
        planet_base["magnetic_field"] = [10.0]

    # Note that the order of Keplerian elements _must not_ be altered.
    if (
        effects["PLANET_TIDES_KAULA_ENABLED"]
        or effects["PLANET_TIDES_CTL_ENABLED"]
        or effects["GENERAL_RELATIVITY_ENABLED"]
    ):
        planet_base.update(
            {
                # No units
                "eccentricity": [0.0],
                # rad
                "pericentre_omega": [0.0],
            }
        )

    if effects["PLANET_TIDES_KAULA_ENABLED"]:
        planet_base.update(
            {
                # rad.s
                "spin": [8.093879511357418e-07],
                # rad
                "inclination": [0.3490658503988659],
                # rad
                "longitude_ascending_node": [1.0],
                # rad
                "spin_inclination": [0.34906584951436426],
                # No units
                "radius_of_gyration": [0.33070368308499226],
                "type_and_file": [
                    (
                        "Solid",
                        "examples/data/planet/tides/kaula/leconte2015_steinberger.json",
                    )
                ],
            }
        )

    if effects["PLANET_TIDES_CTL_ENABLED"]:
        planet_base.update(
            {
                # rad.s-1, initial planet spin (only integrated with the "Free" spin mode)
                "ctl_spin": [8.093879511357418e-07],
                # No units (0.254 for Jupiter)
                "ctl_radius_of_gyration_2": [0.254],
                # s, k2 * delta_t of the planet (Leconte et al. 2010: 2 to 3e-2 s for Jupiter)
                "ctl_k2_delta_t": [0.025],
                # "PseudoSynchronous" (spin held at Leconte et al. 2010 Eq. 12) or "Free"
                "ctl_spin_mode": ["PseudoSynchronous"],
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
        # rad.s-1
        "spin": [TWO_PI / (ROTATORS[REFERENCE][ROTATOR]["period_d"] * 86400.0)],
        # seconds (from years)
        "core_envelope_coupling_constant": [SECONDS_IN_YEAR * x for x in [ROTATORS[REFERENCE][ROTATOR]["coupling_yr"]]],
        "footpoint_conductance": [None],  # Do not edit.
        "evolution": [None],  # Do not edit.
        "sigma_bar": [None],  # Do not edit.
        # J: normalisation of the Matt et al. 2015 wind torque (gamma_sun, Eq. 8).
        # Keep this key after the seven above (they are unpacked by position in config.py).
        "wind_torque_prefactor": [TRACK["prefactor"]],
        "kaula_spectrum_file": [None],  # Do not edit.
    }

    if effects["MAGNETIC_EFFECT_ENABLED"]:
        # Ohm-1
        star_base["footpoint_conductance"] = [5.8e4]

    if effects["STAR_EVOLUTION_ENABLED"]:
        star_base["evolution"] = [TRACK["evolution"]]
    else:
        # Set the initial star values that would otherwise be provided by savgol/mesa data if evolution were enabled.
        # Must be non-zero (to avoid NaN).

        # kg (from Msun)
        star_base["mass"] = [SOLAR_MASS * x for x in [1.0]]
        # m (from SOLAR_RADIUS)
        star_base["radius"] = [SOLAR_RADIUS * x for x in [1.0]]
        # No units
        star_base["radiative_moment_of_inertia"] = [1.0]
        star_base["convective_moment_of_inertia"] = [1.0]

    if effects["STAR_TIDES_KAULA_ENABLED"]:
        star_base["kaula_spectrum_file"] = [SPECTRUM_FILE]

    if effects["STAR_TIDES_CTL_ENABLED"]:
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
            "absolute_tolerance": SECONDS_IN_YEAR * 1e6,
            "relative_tolerance": 1.15,
            "incremental_scaling_factor": 0.002,
            "decremental_scaling_factor": 0.0,
        }
    }

    dopri853 = {
        "Dopri853": {
            "step_controller": {
                "relative_tolerance": 1e-10,
                "absolute_tolerance": 1e-10,
                "step_size_factor_min": 0.3333333333333333,
                "step_size_factor_max": 6.0,
                "step_size_error_factor": 0.9,
                "step_size_max": SECONDS_IN_YEAR * 5e5,
                "alpha": 0.125,
                "beta": 0.0,
            },
            "step_size_underflow": None,
            "stiffness_test": "Disabled",
            "max_integration_steps": 100000000,
            # Uncomment the entire solution_output for Dense output
            "solution_output": filter,
        }
    }

    return dopri853


if __name__ == "__main__":
    make_configs(simulator_setup, effect_setup, planet_setup, star_setup, integrator_setup)
