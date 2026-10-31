"""
Generates dictionaries matching the structure of input configs into spiroid.

"""

import itertools


def make_config(planet, star, disk_lifetime, integrator, start_time, final_time):
    """Generates the simulation template required to launch the simulation in Rust."""
    return {
        "resume": False,
        "initial_time": start_time,
        "final_time": final_time,
        "integrator": integrator,
        "universe": {
            "disk_lifetime": disk_lifetime,
            "central_body": star,
            "orbiting_body": planet,
        },
    }


def make_planets(planet_base, effects):
    """Generate all combinations of planets based on specified values of `planet_base` dictionary."""
    planets = []
    combis = [x for x in itertools.product(*planet_base.values())]

    for planet_vals in combis:
        (mass, radius, semi_major_axis, magnetic_field) = planet_vals[:4]
        body = {}
        planet = {
            "mass": mass,
            "radius": radius,
            "semi_major_axis": semi_major_axis,
        }

        if effects["MAGNETIC_EFFECT_ENABLED"]:
            planet["magnetic_field"] = magnetic_field

        # Remaining entries are looked up by name so that setup scripts written before a
        # key existed keep working.
        named = dict(zip(planet_base.keys(), planet_vals))
        planet_kaula = effects["PLANET_TIDES_KAULA_ENABLED"]
        planet_ctl = effects["PLANET_TIDES_CTL_ENABLED"]
        if planet_kaula and planet_ctl:
            raise ValueError(
                "one tidal model per body: PLANET_TIDES_KAULA_ENABLED and PLANET_TIDES_CTL_ENABLED cannot both be enabled"
            )

        if planet_kaula or planet_ctl or effects["GENERAL_RELATIVITY_ENABLED"]:
            eccentricity = planet_vals[4]
            pericentre_omega = planet_vals[5]
            planet.update({"eccentricity": eccentricity, "pericentre_omega": pericentre_omega})

        # Kaula planetary tide: the planet is the deformed body, its love number comes from
        # the solid file.
        if planet_kaula:
            (
                spin,
                inclination,
                longitude_ascending_node,
                spin_inclination,
                radius_of_gyration,
                (particle_type, kaula_solid_file),
            ) = planet_vals[6:]
            planet.update(
                {
                    "inclination": inclination,
                    "spin": spin,
                    "longitude_ascending_node": longitude_ascending_node,
                    "spin_inclination": spin_inclination,
                    "radius_of_gyration_2": radius_of_gyration,
                }
            )

            body["tides"] = {
                "KaulaTides": {"particle_type": {particle_type: {"solid_file": kaula_solid_file}}}
            }

        # Constant time lag planetary tide (Hut 1981, Leconte et al. 2010): the planet spin
        # (rad.s-1) and radius_of_gyration_2 are the planet inputs, k2_delta_t (s) and the spin
        # mode ("PseudoSynchronous" or "Free") the model inputs.
        if planet_ctl:
            planet.update(
                {
                    "spin": named["ctl_spin"],
                    "radius_of_gyration_2": named["ctl_radius_of_gyration_2"],
                }
            )
            body["tides"] = {
                "PlanetaryConstantTimeLag": {
                    "k2_delta_t": named["ctl_k2_delta_t"],
                    "spin": named["ctl_spin_mode"],
                }
            }

        if not effects["WIND_ENABLED"]:
            body["wind"] = "Disabled"

        body["kind"] = {"Planet": planet}
        planets.append(body)

    return planets


def make_stars(star_base, effects):
    """Generate all combinations of stars based on specified values of `star_base` dictionary."""
    stars = []
    combis = [x for x in itertools.product(*star_base.values())]
    for star_vals in combis:
        (
            mass,
            radius,
            spin,
            core_envelope_coupling_constant,
            footpoint_conductance,
            evolution,
            sigma_bar,
        ) = star_vals[:7]
        # Remaining entries are looked up by name so that setup scripts written before a
        # key existed keep working (the corresponding Rust field then takes its default).
        named = dict(zip(star_base.keys(), star_vals))

        body = {}
        star = {
            "spin": spin,
            "core_envelope_coupling_constant": core_envelope_coupling_constant,
            "evolution": "Disabled",
        }
        # Normalisation of the Matt et al. 2015 wind torque (J): required when the wind is
        # enabled (spiroid refuses to run without it).
        if effects["WIND_ENABLED"]:
            if named.get("wind_torque_prefactor") is None:
                raise KeyError(
                    "star_setup must provide 'wind_torque_prefactor' (J) when WIND_ENABLED"
                )
            star["wind_torque_prefactor"] = named["wind_torque_prefactor"]

        star_kaula = effects["STAR_TIDES_KAULA_ENABLED"]
        star_ctl = effects["STAR_TIDES_CTL_ENABLED"]
        if star_kaula and star_ctl:
            raise ValueError(
                "one tidal model per body: STAR_TIDES_KAULA_ENABLED and STAR_TIDES_CTL_ENABLED cannot both be enabled"
            )

        # Kaula stellar tide: the star is the deformed body, its love number spectrum
        # (examples/data/star/tides/kaula/, named by age, alpha, beta and rotation period)
        # is the star input `kaula_spectrum_file`.
        if star_kaula:
            if named.get("kaula_spectrum_file") is None:
                raise KeyError(
                    "star_setup must provide 'kaula_spectrum_file' when STAR_TIDES_KAULA_ENABLED"
                )
            body["tides"] = {
                "KaulaTides": {
                    "particle_type": {
                        "StarConvectiveEnvelope": {"spectrum_file": named["kaula_spectrum_file"]}
                    }
                }
            }

        # Constant time lag stellar tide (Benbakoura et al. 2019): equilibrium tide with
        # sigma_bar_star plus the frequency-averaged dynamical tide.
        if star_ctl:
            body["tides"] = {
                "StellarConstantTimeLag": {
                    "equilibrium": {"SigmaBarStar": 1e-06},
                    "inertial": "FrequencyAveraged",
                }
            }
        if effects["MAGNETIC_EFFECT_ENABLED"]:
            body["magnetism"] = {"Wind": {"footpoint_conductance": footpoint_conductance}}

        if effects["STAR_EVOLUTION_ENABLED"]:
            star["evolution"] = evolution
        else:
            star["mass"] = mass
            star["radius"] = radius
            star["radiative_moment_of_inertia"] = named["radiative_moment_of_inertia"]
            star["convective_moment_of_inertia"] = named["convective_moment_of_inertia"]

        if not effects["WIND_ENABLED"]:
            body["wind"] = "Disabled"

        if effects.get("GENERAL_RELATIVITY_ENABLED"):
            body["general_relativity"] = "Enabled"

        body["kind"] = {"Star": star}
        stars.append(body)

    return stars


def generate_all_configs(
    start_time, final_time, disk_lifetime, planet_base, star_base, effects, integrator
):
    """Generates a simulation configuration file for each combination of planets and stars."""
    planets = make_planets(planet_base, effects)
    stars = make_stars(star_base, effects)

    # Generate a simulation input config for all combinations
    # of the star and planet values.
    return (
        make_config(planet, star, disk_lifetime, integrator, start_time, final_time)
        for (planet, star) in itertools.product(planets, stars)
    )


def generate_all_effect_combinations(input_dict):
    """Generate all possible combinations of enabled effects for the simulations."""
    import itertools

    tags = {
        "MAGNETIC_EFFECT_ENABLED": "magnetism",
        "STAR_EVOLUTION_ENABLED": "star_evolution",
        "STAR_TIDES_KAULA_ENABLED": "star_kaula_tides",
        "STAR_TIDES_CTL_ENABLED": "star_ctl_tides",
        "PLANET_TIDES_KAULA_ENABLED": "planet_kaula_tides",
        "PLANET_TIDES_CTL_ENABLED": "planet_ctl_tides",
        "WIND_ENABLED": "wind",
        "GENERAL_RELATIVITY_ENABLED": "general_relativity",
    }
    # Get keys and values from the input dictionary
    keys = input_dict.keys()
    values = input_dict.values()
    # Generate all combinations of enabled effects.
    combinations = [x for x in itertools.product(*values)]
    # Create a list of tuples (dictionary, label) where label is the enabled effects.
    result = []
    for combo in combinations:
        # Create the dictionary from the combination
        combo_dict = dict(zip(keys, combo))
        # Concatenate labels from keys with True values (enabled effects) for the filename.
        true_keys = "-".join(tags[key] for key, value in combo_dict.items() if value)
        if true_keys == "":
            true_keys = "no_effects"
        # Append the tuple (dictionary, concatenated string) to the result
        result.append((combo_dict, true_keys))

    return result
