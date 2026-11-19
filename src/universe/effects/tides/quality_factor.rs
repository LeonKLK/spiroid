use crate::constants::GRAVITATIONAL;
use crate::universe::particles::{Planet, Star};

pub(crate) mod zahn;
pub use zahn::Zahn;

pub(crate) mod frequency_averaged;

// Stellar tides expressed with a modified tidal quality factor Q' (Benbakoura et al. 2019):
// the star is the tidally deformed body, on a circular coplanar orbit. These models only
// produce a torque on the star, from which physics.rs derives the stellar spin and the
// semi-major axis (angular momentum conservation); they carry no eccentricity equation.
// - equilibrium tide with the Zahn prescription (`zahn`),
// - frequency-averaged dynamical tide (`frequency_averaged`).
// The constant time lag equilibrium tide is not here: it is `ConstantTimeLag`, valid at any
// eccentricity.

// This is a re-write of Eq. 3 and 19 from Benbakoura et al. 2019
// without the factors that are in the function semi_major_axis_13_div_2_derivative in physics.rs
// The a^-6 is here to compensate the a^6 in physics.rs
// The dissipations of several mechanisms add up (Benbakoura et al. 2019 Eq. 20,
// 1/Q' = 1/Q'_eq + 1/Q'_dyn), so their torques add up too.
pub(crate) fn stellar_tidal_torque(tidal_quality: f64, star: &Star, planet: &Planet) -> f64 {
    // Smoothing parameter when tidal frequency is 0
    let depth = 1e-8;
    -(9. / 4.)
        * planet.mass.powi(2)
        * GRAVITATIONAL
        * planet.semi_major_axis.powi(-6)
        * tanh!(star.tidal_frequency / depth)
        * star.radius.powi(5)
        / tidal_quality
}

#[cfg(test)]
pub mod tests;

// References:
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
