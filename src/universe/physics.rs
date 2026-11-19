use crate::constants::{GRAVITATIONAL, SPEED_OF_LIGHT};
use crate::universe::{
    Kaula, Particle, ParticleType, Planet, PlanetarySpin, Star, UniverseIntegral,
};
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

    // Star derivatives. The stellar tides are not part of these two functions: the tidal
    // torques on the star are collected below (`stellar_tidal_torque`) and added once at
    // the end, the same way the planetary tides are collected for the planet spin.
    dy.central_body.radiative_zone_angular_momentum =
        star_radiative_zone_angular_momentum_derivative(star);
    dy.central_body.convective_zone_angular_momentum =
        star_convective_zone_angular_momentum_derivative(star, disk_is_dissipated);

    // If the planet does not exist, only the star derivatives are computed.
    // i.e. during the disk lifetime, or after the planet is destroyed.
    if !disk_is_dissipated || planet.is_destroyed() {
        return Ok(());
    }

    // The tidal contributions to the spin of each body are first collected over all the
    // enabled tides (equilibrium and dynamical, any model), and only then added to the
    // integrated quantities at the end of this function:
    // - star: sum of the tidal torques on the star (kg.m2.s-2), added to the derivative of
    //   the angular momentum of its convective zone,
    // - planet: sum of the tidal torques on the planet divided by its moment of inertia
    //   (rad.s-2), added to the derivative of its spin.
    // Floating point additions are not associative: adding the tides one by one to the
    // non-tidal terms instead would change the last digit of the result.
    let mut stellar_tidal_torque = 0.0;
    let mut planet_tidal_spin_derivative = 0.0;

    // Quality factor STELLAR tides (`Zahn`, `FrequencyAveraged`; 0 if disabled): a torque on
    // the convective zone of the star. On a circular orbit the angular momentum taken from
    // the stellar spin goes to the orbit, so the semi major axis derivative is obtained from
    // the same torque (angular momentum conservation), together with the magnetic torque.
    // The constant time lag and kaula stellar tides are not in this torque: they write
    // their own semi major axis and eccentricity terms below.
    let stellar_quality_factor_torque = central_body
        .tides
        .stellar_quality_factor_torque(star, planet);
    stellar_tidal_torque += stellar_quality_factor_torque;
    dy.orbiting_body.semi_major_axis = planet_semi_major_axis_13_div_2_derivative_from_torque(
        planet,
        star,
        star.magnetic_torque + stellar_quality_factor_torque,
    );

    // Immutable borrow of kaula properties if kaula planet tides enabled.
    // This block is the PLANETARY tide only (`orbiting_body.tides`): the planet is the
    // deformed body and all `planet_*_derivative` functions below evolve the planet's
    // spin and orbit. A Kaula STELLAR tide (`central_body.tides`) has its own block
    // below and is independent of this one.
    // TODO (stellar tide, non-coplanar): the inclination / node / spin-axis derivatives
    // and the inclination part of the pericentre derivative here are planet-only; see the
    // TODO above `Kaula::summation_of_longitudinal_modes_triple_common`. Coplanar is fine
    // for now.
    if let Some(kaula) = orbiting_body.tides.kaula_get() {
        // Sum the semi major axis derivative to account for both CTL star tide (if enabled) and Kaula planet tide.
        dy.orbiting_body.semi_major_axis +=
            kaula_planet_semi_major_axis_13_div_2_derivative(planet, star, kaula);

        planet_tidal_spin_derivative += planet_spin_derivative(planet, star, kaula);
        dy.orbiting_body.eccentricity = planet_eccentricity_derivative(planet, star, kaula);
        dy.orbiting_body.inclination = planet_inclination_derivative(planet, star, kaula);
        dy.orbiting_body.longitude_ascending_node =
            planet_longitude_ascending_node_derivative(planet, star, kaula);
        dy.orbiting_body.pericentre_omega =
            planet_argument_pericentre_derivative(planet, star, kaula);
        dy.orbiting_body.spin_inclination =
            planet_spin_axis_inclination_derivative(planet, star, kaula);
    }

    // Constant time lag PLANETARY tide (`orbiting_body.tides`, Hut 1981, Leconte et al. 2010):
    // the planet is the deformed body. Coplanar: evolves the semi-major axis, the
    // eccentricity and the planet spin. Angular momentum is exchanged between the planet spin
    // and the orbit only, so nothing is added to the star.
    // With the spin held pseudo-synchronous the torque is not applied: the spin is set from
    // the orbit in `Universe::update_star_planet` instead.
    // Alternative to the kaula planetary tide above (refused together in `Particle::initialise`).
    if let Some(ctl) = orbiting_body.tides.constant_time_lag_get() {
        dy.orbiting_body.semi_major_axis +=
            ctl.semi_major_axis_13_div_2_derivative(planet, star, planet);
        dy.orbiting_body.eccentricity += ctl.eccentricity_squared_derivative(planet, star, planet);
        if ctl.spin == PlanetarySpin::Free {
            planet_tidal_spin_derivative +=
                ctl.tidal_torque(planet, star, planet) / planet.moment_of_inertia;
        }
    }

    // Constant time lag STELLAR tide (`central_body.tides`): the same equations with the
    // star as the deformed body and the planet as the perturber. Evolves the semi-major axis,
    // the eccentricity and the star spin (through the torque on its convective zone).
    if let Some(ctl) = central_body.tides.constant_time_lag_get() {
        dy.orbiting_body.semi_major_axis +=
            ctl.semi_major_axis_13_div_2_derivative(star, planet, planet);
        dy.orbiting_body.eccentricity += ctl.eccentricity_squared_derivative(star, planet, planet);
        stellar_tidal_torque += ctl.tidal_torque(star, planet, planet);
    }

    // Kaula STELLAR tide (`central_body.tides`): the star is the deformed body.
    // Coplanar case only: contributes to the semi major axis, eccentricity and argument
    // of pericentre of the orbit, and its torque to the angular momentum of the star's
    // convective zone.
    // Summed with the planetary tide contributions above (if enabled).
    if let Some(kaula) = central_body.tides.kaula_get() {
        dy.orbiting_body.semi_major_axis +=
            kaula_star_semi_major_axis_13_div_2_derivative(planet, star, kaula);
        dy.orbiting_body.eccentricity += kaula_star_eccentricity_derivative(planet, star, kaula);
        dy.orbiting_body.pericentre_omega +=
            kaula_star_argument_pericentre_derivative(planet, star, kaula);
        stellar_tidal_torque += kaula.tidal_torque_on_deformed_body(star, planet, planet);
    }

    // Add the collected tidal contributions to the spin of each body.
    // Star: the angular momentum of the convective zone is integrated, so the total torque
    // is added as it is. Planet: the spin is integrated (constant moment of inertia).
    dy.central_body.convective_zone_angular_momentum += stellar_tidal_torque;
    dy.orbiting_body.spin += planet_tidal_spin_derivative;

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

// Rate of change in the angular momentum in the convective zone, without the tides.
// Includes additional wind torque which is applicable
// during the post main sequence of the star's evolution.
// The stellar tides are collected and added in `force` (`stellar_tidal_torque`).

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
// `torque_on_star` is the torque applied on the star (kg.m2.s-2, including its a^-6
// dependency) that is balanced by the orbit: magnetism and the quality factor stellar tides.
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
mod tests_mode_ardestani;
#[cfg(test)]
mod tests_mode_fromfile;

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
