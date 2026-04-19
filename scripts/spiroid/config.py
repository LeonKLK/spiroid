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

        if effects["PLANET_TIDES_ENABLED"]:
            (
                spin,
                eccentricity,
                inclination,
                longitude_ascending_node,
                pericentre_omega,
                spin_inclination,
                radius_of_gyration,
                (particle_type, kaula_solid_file),
            ) = planet_vals[4:]
            planet.update(
                {
                    "inclination": inclination,
                    "eccentricity": eccentricity,
                    "spin": spin,
                    "longitude_ascending_node": longitude_ascending_node,
                    "pericentre_omega": pericentre_omega,
                    "spin_inclination": spin_inclination,
                    "radius_of_gyration_2": radius_of_gyration,
                }
            )

            body["tides"] = {
                "KaulaTides": {"particle_type": {particle_type: {"solid_file": kaula_solid_file}}}
            }

        elif effects["STAR_TIDES_ENABLED"] and effects.get("STAR_TIDES_MODEL") == "KaulaTides":
            (
                eccentricity,
                inclination,
                longitude_ascending_node,
                pericentre_omega,
                spin_inclination,
            ) = planet_vals[4:]
            planet.update(
                {
                    "eccentricity": eccentricity,
                    "inclination": inclination,
                    "longitude_ascending_node": longitude_ascending_node,
                    "pericentre_omega": pericentre_omega,
                    "spin_inclination": spin_inclination,
                }
            )

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

        body = {}
        star = {
            "spin": spin,
            "core_envelope_coupling_constant": core_envelope_coupling_constant,
            "evolution": "Disabled",
        }

        if effects["STAR_TIDES_ENABLED"]:
            model = effects.get("STAR_TIDES_MODEL", "ConstantTimeLag")
            if model == "KaulaTides":
                stellar_convective_file = star_vals[7]
                stellar_spectrum_spin_rate = star_vals[8]
                body["tides"] = {
                    "KaulaTides": {
                        "particle_type": {
                            "StellarConvective": {
                                "stellar_convective_file": stellar_convective_file,
                                "stellar_spectrum_spin_rate": stellar_spectrum_spin_rate,
                            }
                        },
                        "atmosphere_model": "Disabled",
                    }
                }
            else:
                body["tides"] = {
                    "ConstantTimeLag": {
                        "equilibrium": {"SigmaBarStar": sigma_bar},
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
            star["radiative_moment_of_inertia"] = star_vals[7]
            star["convective_moment_of_inertia"] = star_vals[8]

        if not effects["WIND_ENABLED"]:
            body["wind"] = "Disabled"

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
        "STAR_TIDES_ENABLED": "star_tides",
        "PLANET_TIDES_ENABLED": "planet_kaula_tides",
        "WIND_ENABLED": "wind",
        "GR_ENABLED": "gr",
    }

    # Separate list values (varying) from scalar values (constants passed through as-is).
    varying = {k: v for k, v in input_dict.items() if isinstance(v, list)}
    constants = {k: v for k, v in input_dict.items() if not isinstance(v, list)}

    keys = list(varying.keys())
    values = list(varying.values())
    combinations = list(itertools.product(*values))

    result = []
    for combo in combinations:
        combo_dict = {**dict(zip(keys, combo)), **constants}
        # Build label from boolean True values; use model name for tides if Kaula.
        label_parts = []
        for key, value in combo_dict.items():
            if key not in tags:
                continue
            if key == "STAR_TIDES_ENABLED" and value:
                model = combo_dict.get("STAR_TIDES_MODEL", "ConstantTimeLag")
                label_parts.append(f"star_{model.lower()}_tides")
            elif value is True:
                label_parts.append(tags[key])
        true_keys = "-".join(label_parts) if label_parts else "no_effects"
        result.append((combo_dict, true_keys))

    return result
