use crate::constants::GRAVITATIONAL;
use crate::universe::particles::{Planet, Star};
use serde::{Deserialize, Serialize};

pub(crate) mod stellar_equilibrium;
pub use stellar_equilibrium::Equilibrium;

pub(crate) mod stellar_inertial;
pub use stellar_inertial::Inertial;

// Stellar constant time lag tide: the star is the tidally deformed body (Benbakoura et al.
// 2019). Produces the tidal torque on the star, from which physics.rs derives the stellar
// spin and the semi-major axis; circular coplanar orbit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub struct StellarConstantTimeLag {
    pub equilibrium: Equilibrium,
    pub inertial: Inertial,
}

impl StellarConstantTimeLag {
    // This is a re-write of Eq. 3 and 19 from Benbakoura et al. 2019
    // without the factors that are in the function semi_major_axis_13_div_2_derivative in physics.rs
    // The a^-6 is here to compensate the a^6 in physics.rs
    pub fn stellar_tidal_torque(&self, star: &Star, planet: &Planet) -> f64 {
        let total_dissipation = 1. / self.equilibrium.tidal_quality(star, planet)
            + 1. / self.inertial.tidal_quality(star, planet);

        let tidal_quality = 1. / total_dissipation;
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
}

// Planetary constant time lag tide: the planet is the tidally deformed body.
//
// Hut (1981) weak-friction model in the form of Leconte et al. (2010), Eqs. 2 to 12: the
// tidal bulge of the planet lags the perturber by a constant time delta_t, so the phase lag
// of every tidal component is proportional to its frequency and the secular equations are
// closed-form polynomials in the eccentricity, valid at any eccentricity. Coplanar case only
// (planet obliquity x_p = cos(eps_p) = 1): the semi-major axis, the eccentricity and the
// planet spin are evolved; angular momentum is exchanged between the planet spin and the
// orbit only, the star is not affected. Unlike `StellarConstantTimeLag`, which produces a
// torque on the star, this model produces the orbit and planet spin derivatives directly
// (see physics.rs).

/// How the planetary constant time lag tide treats the planet spin.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum PlanetarySpin {
    /// The spin is held at its pseudo-synchronous value `omega_eq = N(e) / Omega(e) n`
    /// (Leconte et al. 2010 Eq. 12 with obliquity 0), where the spin derivative vanishes.
    /// The spin relaxes to `omega_eq` on a timescale `C_p n^2 / K_p` (about 1e5 yr for a
    /// Jupiter at 0.05 AU, 100 yr at 0.015 AU), which makes the integration of the spin stiff
    /// over Gyr while the orbit does not feel the difference once the spin is locked.
    #[default]
    PseudoSynchronous,
    /// The spin is integrated from its input value (Eq. 9). Needs the planet input
    /// `radius_of_gyration_2` for its moment of inertia.
    Free,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields, default)]
pub struct PlanetaryConstantTimeLag {
    /// Product of the planet's second order potential love number and its constant time lag,
    /// `k2_p * delta_t_p` (s). Leconte et al. (2010): 2 to 3e-2 s for Jupiter; 213 s for the
    /// Earth (Neron de Surgy & Laskar 1997). Must be positive.
    pub k2_delta_t: f64,
    pub spin: PlanetarySpin,
}

impl PlanetaryConstantTimeLag {
    /// Tidal prefactor `K_p` (W), Leconte et al. (2010) Eq. 5:
    /// ```text
    /// K_p = 3/2 k2_p dt_p (G M_p^2 / R_p) (M_star / M_p)^2 (R_p / a)^6 n^2
    ///     = 3/2 k2_p dt_p G M_star^2 R_p^5 n^2 / a^6
    /// ```
    pub(crate) fn prefactor(&self, planet: &Planet, star: &Star) -> f64 {
        1.5 * self.k2_delta_t
            * GRAVITATIONAL
            * star.mass.powi(2)
            * planet.radius.powi(5)
            * planet.mean_motion.powi(2)
            / planet.semi_major_axis.powi(6)
    }

    // Eccentricity functions of Hut (1981), numbered as in Leconte et al. (2010).
    // All equal 1 for a circular orbit.

    /// `N(e)`, Eq. 3.
    pub(crate) fn big_n(e: f64) -> f64 {
        let e2 = e * e;
        (1. + 15. / 2. * e2 + 45. / 8. * e2 * e2 + 5. / 16. * e2 * e2 * e2) / (1. - e2).powi(6)
    }

    /// `N_a(e)`, Eq. 4.
    pub(crate) fn big_n_a(e: f64) -> f64 {
        let e2 = e * e;
        (1. + 31. / 2. * e2
            + 255. / 8. * e2 * e2
            + 185. / 16. * e2 * e2 * e2
            + 25. / 64. * e2 * e2 * e2 * e2)
            / (1. - e2).powf(7.5)
    }

    /// `Omega_e(e)`, Eq. 7.
    pub(crate) fn big_omega_e(e: f64) -> f64 {
        let e2 = e * e;
        (1. + 3. / 2. * e2 + 1. / 8. * e2 * e2) / (1. - e2).powi(5)
    }

    /// `N_e(e)`, Eq. 8.
    pub(crate) fn big_n_e(e: f64) -> f64 {
        let e2 = e * e;
        (1. + 15. / 4. * e2 + 15. / 8. * e2 * e2 + 5. / 64. * e2 * e2 * e2) / (1. - e2).powf(6.5)
    }

    /// `Omega(e)`, Eq. 11.
    pub(crate) fn big_omega(e: f64) -> f64 {
        let e2 = e * e;
        (1. + 3. * e2 + 3. / 8. * e2 * e2) / (1. - e2).powf(4.5)
    }

    /// Pseudo-synchronous spin of the planet, Eq. 12 with obliquity 0:
    /// `omega_eq = N(e) / Omega(e) n`. The spin derivative vanishes there; equals `n` for `e = 0`.
    pub(crate) fn pseudo_synchronous_spin(planet: &Planet) -> f64 {
        Self::big_n(planet.eccentricity) / Self::big_omega(planet.eccentricity) * planet.mean_motion
    }

    /// Tidal torque exerted on the planet by the star (kg.m2.s-2), Eq. 9 with `x_p = 1`:
    /// ```text
    /// C_p d(omega_p)/dt = -(K_p / n) [2 Omega(e) omega_p / n - 2 N(e)]
    /// ```
    /// For `e = 0` this is Hut's `-3 k2 dt G M_star^2 R_p^5 (omega_p - n) / a^6`.
    pub(crate) fn planetary_tidal_torque(&self, planet: &Planet, star: &Star) -> f64 {
        let e = planet.eccentricity;
        -self.prefactor(planet, star) / planet.mean_motion
            * (2. * Self::big_omega(e) * planet.spin / planet.mean_motion - 2. * Self::big_n(e))
    }

    /// Derivative of `a^(13/2)` (the integrated variable) from the planetary tide.
    /// Leconte et al. (2010) Eq. 2, planetary term:
    /// ```text
    /// (1/a) da/dt = 4 a / (G M_star M_p) K_p [N(e) omega_p / n - N_a(e)]
    /// ```
    pub(crate) fn semi_major_axis_13_div_2_derivative(&self, planet: &Planet, star: &Star) -> f64 {
        let e = planet.eccentricity;
        let relative_rate = 4. * planet.semi_major_axis / (GRAVITATIONAL * star.mass * planet.mass)
            * self.prefactor(planet, star)
            * (Self::big_n(e) * planet.spin / planet.mean_motion - Self::big_n_a(e));

        6.5 * planet.semi_major_axis.powf(6.5) * relative_rate
    }

    /// Derivative of `e^2` (the integrated variable) from the planetary tide.
    /// Leconte et al. (2010) Eq. 6, planetary term:
    /// ```text
    /// (1/e) de/dt = 11 a / (G M_star M_p) K_p [Omega_e(e) omega_p / n - 18/11 N_e(e)]
    /// ```
    /// so `d(e^2)/dt = 2 e^2 (1/e) de/dt`, which is 0 for a circular orbit.
    pub(crate) fn eccentricity_squared_derivative(&self, planet: &Planet, star: &Star) -> f64 {
        let e = planet.eccentricity;
        if e == 0. {
            return 0.;
        }
        let relative_rate = 11. * planet.semi_major_axis
            / (GRAVITATIONAL * star.mass * planet.mass)
            * self.prefactor(planet, star)
            * (Self::big_omega_e(e) * planet.spin / planet.mean_motion
                - 18. / 11. * Self::big_n_e(e));

        2. * e * e * relative_rate
    }

    /// Derivative of the planet spin (rad.s-2): the torque divided by the planet's moment
    /// of inertia (constant). Zero when the spin is held pseudo-synchronous.
    pub(crate) fn spin_derivative(&self, planet: &Planet, star: &Star) -> f64 {
        match self.spin {
            PlanetarySpin::PseudoSynchronous => 0.,
            PlanetarySpin::Free => {
                self.planetary_tidal_torque(planet, star) / planet.moment_of_inertia
            }
        }
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Hut 1981, https://ui.adsabs.harvard.edu/abs/1981A%26A....99..126H/abstract
// Leconte et al. 2010, https://doi.org/10.1051/0004-6361/201014337
