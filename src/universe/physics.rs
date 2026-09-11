use crate::constants::{GRAVITATIONAL, SPEED_OF_LIGHT};
use crate::universe::effects::tides::TidalModel;
use crate::universe::{Kaula, Particle, ParticleType, Planet, Star, UniverseIntegral};
use anyhow::{Result, bail};

pub(crate) fn force(
    central_body: &Particle,
    orbiting_body: &Particle,
    disk_is_dissipated: bool,
    dy: &mut UniverseIntegral,
) -> Result<()> {
    dy.zero();

    let ParticleType::Star(star) = &central_body.kind else {
        todo!();
    };

    let ParticleType::Planet(planet) = &orbiting_body.kind else {
        todo!();
    };

    // Star derivatives
    dy.central_body.radiative_zone_angular_momentum =
        star_radiative_zone_angular_momentum_derivative(star);
    dy.central_body.convective_zone_angular_momentum =
        star_convective_zone_angular_momentum_derivative(star, disk_is_dissipated);

    // If the planet does not exist, only the star derivatives are computed.
    // i.e. during the disk lifetime, or after the planet is destroyed.
    if !disk_is_dissipated || planet.is_destroyed() {
        return Ok(());
    }

    // Semi major axis derivative from the torques applied on the star (angular momentum
    // conservation): magnetism and the constant time lag stellar tide (0 if disabled).
    // For a kaula stellar tide the torque only evolves the star's spin (see
    // `star_convective_zone_angular_momentum_derivative`); its effect on the orbit is
    // taken from the kaula summations below instead, which also carry the eccentricity
    // terms, so the tidal torque is excluded here to not count it twice.
    let star_kaula = central_body.tides.kaula_get();
    dy.orbiting_body.semi_major_axis = match star_kaula {
        None => planet_semi_major_axis_13_div_2_derivative(planet, star),
        Some(_) => planet_semi_major_axis_13_div_2_derivative_from_torque(
            planet,
            star,
            star.magnetic_torque,
        ),
    };

    // Immutable borrow of kaula properties if kaula planet tides enabled.
    // This block is the PLANETARY tide only (`orbiting_body.tides`): the planet is the
    // deformed body and all `planet_*_derivative` functions below evolve the planet's
    // spin and orbit. A Kaula STELLAR tide (`central_body.tides`) enters through
    // `star.tidal_torque_convective` instead (see `TidalModel::tidal_torque`), and is
    // independent of this block.
    // TODO (stellar tide, non-coplanar): the inclination / node / spin-axis derivatives
    // and the inclination part of the pericentre derivative here are planet-only; see the
    // TODO above `Kaula::summation_of_longitudinal_modes_triple_common`. Coplanar is fine
    // for now.
    if let TidalModel::KaulaTides(ref kaula) = orbiting_body.tides {
        // Sum the semi major axis derivative to account for both CTL star tide (if enabled) and Kaula planet tide.
        dy.orbiting_body.semi_major_axis +=
            kaula_planet_semi_major_axis_13_div_2_derivative(planet, star, kaula);

        dy.orbiting_body.spin = planet_spin_derivative(planet, star, kaula);
        dy.orbiting_body.eccentricity = planet_eccentricity_derivative(planet, star, kaula);
        dy.orbiting_body.inclination = planet_inclination_derivative(planet, star, kaula);
        dy.orbiting_body.longitude_ascending_node =
            planet_longitude_ascending_node_derivative(planet, star, kaula);
        dy.orbiting_body.pericentre_omega =
            planet_argument_pericentre_derivative(planet, star, kaula);
        dy.orbiting_body.spin_inclination =
            planet_spin_axis_inclination_derivative(planet, star, kaula);
    }

    // Kaula STELLAR tide (`central_body.tides`): the star is the deformed body.
    // Coplanar case only: contributes to the semi major axis, eccentricity and argument
    // of pericentre of the orbit; the star's spin is evolved through
    // `star.tidal_torque_convective`.
    // Summed with the planetary tide contributions above (if enabled).
    if let Some(kaula) = star_kaula {
        dy.orbiting_body.semi_major_axis +=
            kaula_star_semi_major_axis_13_div_2_derivative(planet, star, kaula);
        dy.orbiting_body.eccentricity += kaula_star_eccentricity_derivative(planet, star, kaula);
        dy.orbiting_body.pericentre_omega +=
            kaula_star_argument_pericentre_derivative(planet, star, kaula);
    }

    // General Relativity 1PN apsidal precession.
    if central_body.general_relativity.is_enabled() {
        dy.orbiting_body.pericentre_omega += general_relativity_pericentre_precession_rate(planet);
    }

    // Check the derivatives for numerical errors.
    if dy.denormal_check() {
        let msg = format!("{:?}, {:?}, dy: {:?}", &star, &planet, &dy);
        eprintln!("{}", &msg);
        bail!(
            "error in computation of derivatives: Houston, we have a NaN...infinity, and beyond! {msg}"
        );
    }

    Ok(())
}

// Rate of change in the angular momentum in the convective zone.
// Includes additional wind torque which is applicable
// during the post main sequence of the star's evolution.

// Ahuir et al. 2021, Eq. 2.
fn star_convective_zone_angular_momentum_derivative(star: &Star, disk_is_dissipated: bool) -> f64 {
    // If the disk has not dissipated, the spin of the star has not evolved.
    if !disk_is_dissipated {
        return star.spin * star.convective_moment_of_inertia_derivative;
    }
    star.angular_momentum_redistribution / star.core_envelope_coupling_constant
        - star.mass_transfer_envelope_to_core_torque
        + star.wind_torque
        // evolved_wind_torque should be zero if not in the post main sequence.
        + star.evolved_wind_torque
        + star.magnetic_torque
        + star.tidal_torque_convective
}

// Rate of change in the angular momentum in the radiative zone.
// Benbakoura et al. 2019, Eq. 2
fn star_radiative_zone_angular_momentum_derivative(star: &Star) -> f64 {
    -star.angular_momentum_redistribution / star.core_envelope_coupling_constant
        + star.mass_transfer_envelope_to_core_torque
}

// This loosely comes from Eq. 1 from Ahuir et al. 2021 (for the sum of tidal and magnetic components)
// Also Benbakoura et al. 2019, Eq. 3.
// It is the derivative of semi major axis (a) to the power 13/2
// This is obtained by moving the 1/a^6 dependency of the tidal torque to the left of Eq. 3, alongside the a^(1/2)
// this means that what we call here the tidal torque is not exactly the tidal torque, but the tidal torque * a^6
// or tidal torque without the semi-major axis dependency
// The last line corresponds to the change in semi-major axis from the mass lost in the evolved phases of evolution.
// evolved_change_semi_major_axis is da/dt so is multiplied by 13/2 a^{11/2} to represent the derivative of a^{13/2}
pub(crate) fn planet_semi_major_axis_13_div_2_derivative(planet: &Planet, star: &Star) -> f64 {
    planet_semi_major_axis_13_div_2_derivative_from_torque(
        planet,
        star,
        star.magnetic_torque + star.tidal_torque_convective,
    )
}

// Same as `planet_semi_major_axis_13_div_2_derivative` with an explicit torque applied
// on the star (kg.m2.s-2, including its a^-6 dependency).
pub(crate) fn planet_semi_major_axis_13_div_2_derivative_from_torque(
    planet: &Planet,
    star: &Star,
    torque_on_star: f64,
) -> f64 {
    -13. * sqrt!((star.mass + planet.mass) / GRAVITATIONAL)
        * (1. / (star.mass * planet.mass))
        * planet.semi_major_axis.powi(6)
        * torque_on_star
        + 13. / 2. * planet.semi_major_axis.powf(11. / 2.) * star.evolved_change_semi_major_axis
}

// Semi-major axis derivative from the kaula STELLAR tide (star = deformed body).
// `kaula_planet_semi_major_axis_13_div_2_derivative` with the roles of star and planet
// swapped: radius of the deformed body (star), mass ratio perturber / deformed.
// For e = 0 and i = 0 this equals the torque route
// (`planet_semi_major_axis_13_div_2_derivative_from_torque` with
// `Kaula::tidal_torque_on_deformed_body`), since the spin and semi-major axis
// summations coincide when only (m, p, q) = (2, 0, 0) contributes.
pub(crate) fn kaula_star_semi_major_axis_13_div_2_derivative(
    planet: &Planet,
    star: &Star,
    kaula: &Kaula,
) -> f64 {
    -13. * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
        * (planet.mass / star.mass)
        * star.radius.powi(5)
        * kaula.summation_of_longitudinal_modes_semi_major_axis()
}

// Eccentricity derivative from the kaula STELLAR tide (star = deformed body).
// `planet_eccentricity_derivative` with the roles of star and planet swapped.
fn kaula_star_eccentricity_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if planet.eccentricity == 0. {
        0.
    } else {
        -2.0 * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
            * (star.radius.powi(5) / planet.semi_major_axis.powf(6.5))
            * (planet.mass / star.mass)
            * planet.semi_minor_axis_ratio
            * kaula.summation_of_longitudinal_modes_eccentricity()
    }
}

// Argument of pericentre derivative from the kaula STELLAR tide (star = deformed body).
// `planet_argument_pericentre_derivative` with the roles of star and planet swapped:
// prefactor G M_planet^2 R_star^5 / a^6 (the torque prefactor of the deformed star).
// Coplanar case only: the eccentricity part (Boue & Efroimksy (2019) Eq. 120, Revol et
// al. (2023) Eq A.11) only reads orbital quantities and is role-agnostic. The inclination
// part needs the stellar spin axis and is not implemented (see the TODO above
// `Kaula::summation_of_longitudinal_modes_triple_common`); it is zero in the coplanar case.
fn kaula_star_argument_pericentre_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if planet.eccentricity == 0. {
        0.
    } else {
        ((GRAVITATIONAL * planet.mass.powi(2) * star.radius.powi(5))
            / planet.semi_major_axis.powi(6))
            * kaula.summation_of_longitudinal_modes_pericentre_eccentricity(planet)
    }
}

// Semi-major axis derivative.
// Boue & Efroimksy (2019) Eq. 116 and Revol et al. (2023) Eq A.1
pub(crate) fn kaula_planet_semi_major_axis_13_div_2_derivative(
    planet: &Planet,
    star: &Star,
    kaula: &Kaula,
) -> f64 {
    -13. * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
        * (star.mass / planet.mass)
        * planet.radius.powi(5)
        * kaula.summation_of_longitudinal_modes_semi_major_axis()
}

// Spin derivative.
// Boue & Efroimksy (2019) Eq. 123 and Revol et al. (2023) Eq A.3
fn planet_spin_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    let planet_tidal_torque = (GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5))
        / (planet.semi_major_axis.powi(6));

    (planet_tidal_torque / planet.moment_of_inertia) * kaula.summation_of_longitudinal_modes_spin()
}

// Eccentricity derivative.
// Boue & Efroimksy (2019) Eq. 117 and Revol et al. (2023) Eq A.3
fn planet_eccentricity_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if planet.eccentricity == 0. {
        0.
    } else {
        -2.0 * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
            * (planet.radius.powi(5) / planet.semi_major_axis.powf(6.5))
            * (star.mass / planet.mass)
            * planet.semi_minor_axis_ratio
            * kaula.summation_of_longitudinal_modes_eccentricity()
    }
}

// Inclination derivative.
// Boue & Efroimksy (2019) Eq. 118 and Revol et al. (2023) Eq A.7
// The inclination refers to the angle between the orbital planet and the planet's equatorial plane (i.e. obliquity)
fn planet_inclination_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if planet.sin_inc == 0.0 {
        0.0
    } else {
        (1. / planet.sin_inc)
            * (star.mass / planet.mass)
            * (planet.radius / planet.semi_major_axis).powi(5)
            * kaula.summation_of_longitudinal_modes_inclination()
    }
}

// Longitude of ascending node derivative.
// Boue & Efroimksy (2019) Eq. 121 and Revol et al. (2023) Eq A.9
fn planet_longitude_ascending_node_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if (planet.inclination == 0.) || (planet.spin_inclination == 0.) {
        0.0
    } else {
        ((GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5))
            / planet.semi_major_axis.powi(6))
            * kaula.summation_of_longitudinal_modes_longitude_ascending_node(planet)
    }
}

// Longitude of pericentre derivative.
// Boue & Efroimksy (2019) Eq. 120 and Revol et al. (2023) Eq A.11
fn planet_argument_pericentre_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    // inclination
    let summation_of_longitudinal_modes_pericentre_inclination =
        if planet.inclination == 0. || planet.spin_inclination == 0. {
            0.
        } else {
            kaula.summation_of_longitudinal_modes_pericentre_inclination(planet)
        };

    // eccentricity
    let summation_of_longitudinal_modes_pericentre_eccentricity = if planet.eccentricity == 0. {
        0.
    } else {
        kaula.summation_of_longitudinal_modes_pericentre_eccentricity(planet)
    };

    ((GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5)) / planet.semi_major_axis.powi(6))
        * (summation_of_longitudinal_modes_pericentre_eccentricity
            + summation_of_longitudinal_modes_pericentre_inclination)
}

// Spin axis inclination derivative.
// Boue & Efroimksy (2019) Eq. 122 and Revol et al. (2023) Eq. A.13
// The spin axis inclination refers to the inclination of the planet's rotational vector
// with respect to the total angular momentum.
fn planet_spin_axis_inclination_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if (planet.inclination == 0.0) || (planet.spin_inclination == 0.0) {
        0.
    } else {
        (GRAVITATIONAL * star.mass.powi(2) * planet.radius.powi(5))
            / (planet.semi_major_axis.powi(6) * planet.moment_of_inertia * planet.spin)
            * kaula.summation_of_longitudinal_modes_spin_axis_inclination(planet)
    }
}

#[cfg(test)]
mod tests;

// General Relativity 1PN apsidal precession rate.
// Einstein (1915)
// Weinberg (1972) Eq. 9.2.5
// dω/dt = 3 n³ a² / (c² (1 - e²))
// Derived from dω/dt = 3 G M n / (c² a (1 - e²)) via Kepler's 3rd law: GM = n² a³
fn general_relativity_pericentre_precession_rate(planet: &Planet) -> f64 {
    3.0 * planet.mean_motion.powi(3) * planet.semi_major_axis.powi(2)
        / (SPEED_OF_LIGHT.powi(2) * (1.0 - planet.eccentricity.powi(2)))
}

// References:
// Ahuir et al. 2021, https://doi.org/10.1051/0004-6361/202040173
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Boué and Efroimsky 2019, https://doi.org/10.1007/s10569-019-9908-2
// Einstein 1915, https://doi.org/10.1002/3527608958.ch7
// Revol et al. 2023, https://doi.org/10.1051/0004-6361/202245790
// Weinberg 1972, Gravitation and Cosmology, Eq. 9.2.5
