use crate::constants::GRAVITATIONAL;
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

    // Star tidal contributions — dispatch on the star's chosen tidal model.
    // Each branch is responsible for including non-tidal contributions (magnetic, evolved mass loss)
    // to preserve floating-point order and keep golden file outputs stable.
    match &central_body.tides {
        TidalModel::Disabled => {
            dy.orbiting_body.semi_major_axis +=
                planet_semi_major_axis_13_div_2_non_tidal(planet, star);
        }
        TidalModel::ConstantTimeLag(_) => {
            dy.orbiting_body.semi_major_axis +=
                planet_semi_major_axis_13_div_2_derivative(planet, star);
        }
        TidalModel::KaulaTides(kaula) => {
            dy.orbiting_body.semi_major_axis +=
                planet_semi_major_axis_13_div_2_non_tidal(planet, star)
                    + kaula_star_semi_major_axis_13_div_2_tidal(planet, star, kaula);
            dy.central_body.convective_zone_angular_momentum +=
                kaula_star_convective_zone_angular_momentum_derivative(planet, star, kaula);
            dy.orbiting_body.eccentricity +=
                kaula_star_eccentricity_derivative(planet, star, kaula);
            dy.orbiting_body.inclination +=
                kaula_star_inclination_derivative(planet, star, kaula);
            dy.orbiting_body.longitude_ascending_node +=
                kaula_star_longitude_ascending_node_derivative(planet, star, kaula);
            dy.orbiting_body.pericentre_omega +=
                kaula_star_argument_pericentre_derivative(planet, star, kaula);
            dy.orbiting_body.spin_inclination +=
                kaula_star_spin_axis_inclination_derivative(planet, star, kaula);
        }
    }

    // Planet tidal contributions — dispatch on the planet's chosen tidal model.
    match &orbiting_body.tides {
        TidalModel::Disabled => {}
        TidalModel::ConstantTimeLag(_) => todo!("planet CTL tides"),
        TidalModel::KaulaTides(kaula) => {
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

// Non-tidal orbital contributions (magnetic torque + evolved mass loss) with no stellar tide.
// Used when the star's tidal model is Disabled.
// Ahuir et al. 2021, Eq. 1 (magnetic component) and mass-loss term.
// evolved_change_semi_major_axis is da/dt, multiplied by 13/2 a^{11/2} to get d(a^{13/2})/dt.
fn planet_semi_major_axis_13_div_2_non_tidal(planet: &Planet, star: &Star) -> f64 {
    -13. * sqrt!((star.mass + planet.mass) / GRAVITATIONAL)
        * (1. / (star.mass * planet.mass))
        * planet.semi_major_axis.powi(6)
        * star.magnetic_torque
        + 13. / 2. * planet.semi_major_axis.powf(11. / 2.) * star.evolved_change_semi_major_axis
}

// CTL stellar tide + non-tidal contributions combined into a single expression.
// Preserving the combined form keeps floating-point order identical to the original.
// Ahuir et al. 2021, Eq. 1; Benbakoura et al. 2019, Eq. 3.
// The a^{-6} dependency of the tidal torque is moved to the left alongside a^{1/2},
// so tidal_torque_convective here excludes the semi-major axis dependency.
pub(crate) fn planet_semi_major_axis_13_div_2_derivative(planet: &Planet, star: &Star) -> f64 {
    -13. * sqrt!((star.mass + planet.mass) / GRAVITATIONAL)
        * (1. / (star.mass * planet.mass))
        * planet.semi_major_axis.powi(6)
        * (star.magnetic_torque + star.tidal_torque_convective)
        + 13. / 2. * planet.semi_major_axis.powf(11. / 2.) * star.evolved_change_semi_major_axis
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
// Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.13
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

// --- Star Kaula tidal functions ---
// Mirror of the planet Kaula functions with roles swapped:
// star is the deformed body, planet is the perturber.

// Semi-major axis tidal contribution from star Kaula tide.
// Boue & Efroimksy (2019) Eq. 116 and Revol et al. (2023) Eq A.1
fn kaula_star_semi_major_axis_13_div_2_tidal(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    -13. * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
        * (planet.mass / star.mass)
        * star.radius.powi(5)
        * kaula.summation_of_longitudinal_modes_semi_major_axis()
}

// Tidal torque on the star's convective zone from Kaula tide.
// Boue & Efroimksy (2019) Eq. 123 and Revol et al. (2023) Eq A.3
// For star we integrate the angular momentum instead of the spin rate
fn kaula_star_convective_zone_angular_momentum_derivative(
    planet: &Planet,
    star: &Star,
    kaula: &Kaula,
) -> f64 {
    let star_tidal_torque = (GRAVITATIONAL * planet.mass.powi(2) * star.radius.powi(5))
        / star.kaula_semi_major_axis.powi(6);
    
    star_tidal_torque * kaula.summation_of_longitudinal_modes_spin()
}

// Eccentricity derivative from star Kaula tide.
// Boue & Efroimksy (2019) Eq. 117 and Revol et al. (2023) Eq A.3
fn kaula_star_eccentricity_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if star.kaula_eccentricity == 0. {
        0.
    } else {
        -2.0 * sqrt!(GRAVITATIONAL * (star.mass + planet.mass))
            * (star.radius.powi(5) / star.kaula_semi_major_axis.powf(6.5))
            * (planet.mass / star.mass)
            * star.kaula_semi_minor_axis_ratio
            * kaula.summation_of_longitudinal_modes_eccentricity()
    }
}

// Inclination derivative from star Kaula tide.
// Boue & Efroimksy (2019) Eq. 118 and Revol et al. (2023) Eq A.7
fn kaula_star_inclination_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if star.kaula_sin_inc == 0.0 {
        0.0
    } else {
        (1. / star.kaula_sin_inc)
            * (planet.mass / star.mass)
            * (star.radius / star.kaula_semi_major_axis).powi(5)
            * kaula.summation_of_longitudinal_modes_inclination()
    }
}

// Longitude of ascending node derivative from star Kaula tide.
// Boue & Efroimksy (2019) Eq. 121 and Revol et al. (2023) Eq A.9
fn kaula_star_longitude_ascending_node_derivative(
    planet: &Planet,
    star: &Star,
    kaula: &Kaula,
) -> f64 {
    if star.kaula_inclination == 0. || star.kaula_spin_inclination == 0. {
        0.0
    } else {
        ((GRAVITATIONAL * planet.mass.powi(2) * star.radius.powi(5))
            / star.kaula_semi_major_axis.powi(6))
            * kaula.summation_of_longitudinal_modes_longitude_ascending_node_stellar(star)
    }
}

// Longitude of pericentre derivative from star Kaula tide.
// Boue & Efroimksy (2019) Eq. 120 and Revol et al. (2023) Eq A.11
fn kaula_star_argument_pericentre_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    let summation_pericentre_inclination =
        if star.kaula_inclination == 0. || star.kaula_spin_inclination == 0. {
            0.
        } else {
            kaula.summation_of_longitudinal_modes_pericentre_inclination_stellar(star)
        };

    let summation_pericentre_eccentricity = if star.kaula_eccentricity == 0. {
        0.
    } else {
        kaula.summation_of_longitudinal_modes_pericentre_eccentricity_stellar(star)
    };

    ((GRAVITATIONAL * planet.mass.powi(2) * star.radius.powi(5))
        / star.kaula_semi_major_axis.powi(6))
        * (summation_pericentre_eccentricity + summation_pericentre_inclination)
}

// Spin axis inclination derivative from star Kaula tide.
// Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.13
fn kaula_star_spin_axis_inclination_derivative(planet: &Planet, star: &Star, kaula: &Kaula) -> f64 {
    if star.kaula_inclination == 0.0 || star.kaula_spin_inclination == 0.0 {
        0.
    } else {
        let moment_of_inertia =
            star.convective_moment_of_inertia + star.radiative_moment_of_inertia;
        (GRAVITATIONAL * planet.mass.powi(2) * star.radius.powi(5))
            / (star.kaula_semi_major_axis.powi(6) * moment_of_inertia * star.spin)
            * kaula.summation_of_longitudinal_modes_spin_axis_inclination_stellar(star)
    }
}

#[cfg(test)]
mod tests;

// References:
// Ahuir et al. 2021, https://doi.org/10.1051/0004-6361/202040173
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Boué and Efroimsky 2019, https://doi.org/10.1007/s10569-019-9908-2
// Revol et al. 2023, https://doi.org/10.1051/0004-6361/202245790
