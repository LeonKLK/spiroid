use crate::constants::{GRAVITATIONAL, SOLAR_MASS, SOLAR_RADIUS};
use crate::universe::particles::{ParticleT, Planet};
use serde::{Deserialize, Serialize};

// Constant time lag equilibrium tide, shared by the star and the planet.
//
// Hut (1981) weak-friction model in the form of Leconte et al. (2010), Eqs. 2 to 12: the
// tidal bulge of the deformed body lags the perturber by a constant time delta_t, so the
// phase lag of every tidal component is proportional to its frequency and the secular
// equations are closed-form polynomials in the eccentricity, valid at any eccentricity.
// Coplanar case only (obliquity x = cos(eps) = 1).
//
// The same equations serve both bodies. Every function takes the roles explicitly:
// - `deformed`: the body that carries the tidal bulge (its radius, its spin),
// - `perturber`: the body that raises the tide (its mass),
// - `orbit`: the body that carries the orbital elements (always the planet: a, n, e).
// Stellar tide: `(star, planet, planet)`. Planetary tide: `(planet, star, planet)`.
//
// The model evolves the semi-major axis, the eccentricity and the spin of the deformed
// body; angular momentum is exchanged between that spin and the orbit only (see physics.rs).

/// How the constant time lag tide treats the spin of the PLANET (ignored on the star, whose
/// spin is always integrated through the angular momentum of its convective zone).
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
pub struct ConstantTimeLag {
    /// Product of the second order potential love number of the deformed body and its
    /// constant time lag, `k2 * delta_t` (s), constant in time. Leconte et al. (2010): 2 to
    /// 3e-2 s for Jupiter; 213 s for the Earth (Neron de Surgy & Laskar 1997).
    /// Give either this or `sigma_bar_star`.
    pub k2_delta_t: f64,
    /// Star only: normalised dissipation of Hansen (2010, 2012), `sigma_bar_star`
    /// (dimensionless, e.g. 7.8e-8), constant in time. The time lag then follows the stellar
    /// radius, see `ConstantTimeLag::k2_delta_t`. Give either this or `k2_delta_t`.
    pub sigma_bar_star: f64,
    /// Planet only, see `PlanetarySpin`.
    pub spin: PlanetarySpin,
}

impl ConstantTimeLag {
    /// `k2 * delta_t` (s) of the deformed body.
    ///
    /// Either the constant input `k2_delta_t`, or, from the normalised dissipation
    /// `sigma_bar_star` (Bolmont & Mathis 2016 Eq. 8, Bolmont et al. 2015 Eq. 28):
    /// ```text
    /// k2 dt = 3/2 sigma_bar_star sigma_0 R^5 / G,    sigma_0 = sqrt(G / (M_sun R_sun^7))
    /// ```
    /// i.e. 2389 s * `sigma_bar_star` * (R / `R_sun`)^5. This is the relation that reproduces
    /// the definition of the dissipation by the orbital decay rate (Hansen 2010 Eq. 5,
    /// Bolmont & Mathis 2016 Eqs. 12 and 13). Benbakoura et al. (2019) Eq. 21 writes it
    /// without the factor 3/2, which is where the extra 3/2 of their Eq. 22 comes from.
    pub(crate) fn k2_delta_t(&self, deformed: &impl ParticleT) -> f64 {
        if self.sigma_bar_star > 0. {
            let normalisation_constant = sqrt!(GRAVITATIONAL / (SOLAR_MASS * SOLAR_RADIUS.powi(7)));
            1.5 * self.sigma_bar_star * normalisation_constant * deformed.radius().powi(5)
                / GRAVITATIONAL
        } else {
            self.k2_delta_t
        }
    }

    /// Tidal prefactor `K` (W) of the deformed body, Leconte et al. (2010) Eq. 5:
    /// ```text
    /// K = 3/2 k2 dt (G M_def^2 / R_def) (M_pert / M_def)^2 (R_def / a)^6 n^2
    ///   = 3/2 k2 dt G M_pert^2 R_def^5 n^2 / a^6
    /// ```
    pub(crate) fn prefactor(
        &self,
        deformed: &impl ParticleT,
        perturber: &impl ParticleT,
        orbit: &Planet,
    ) -> f64 {
        1.5 * self.k2_delta_t(deformed)
            * GRAVITATIONAL
            * perturber.mass().powi(2)
            * deformed.radius().powi(5)
            * orbit.mean_motion.powi(2)
            / orbit.semi_major_axis.powi(6)
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

    /// Tidal torque exerted on the deformed body by the perturber (kg.m2.s-2), Eq. 9 with
    /// `x = 1`:
    /// ```text
    /// C d(omega)/dt = -(K / n) [2 Omega(e) omega / n - 2 N(e)]
    /// ```
    /// For `e = 0` this is Hut's `-3 k2 dt G M_pert^2 R_def^5 (omega - n) / a^6`, which for the
    /// star is the torque of Benbakoura et al. (2019) Eq. 19 with the equilibrium tide
    /// quality factor `Q' = G / (sigma_star R^5 |2 (omega - n)|)`.
    pub(crate) fn tidal_torque(
        &self,
        deformed: &impl ParticleT,
        perturber: &impl ParticleT,
        orbit: &Planet,
    ) -> f64 {
        let e = orbit.eccentricity;
        -self.prefactor(deformed, perturber, orbit) / orbit.mean_motion
            * (2. * Self::big_omega(e) * deformed.spin() / orbit.mean_motion - 2. * Self::big_n(e))
    }

    /// Derivative of `a^(13/2)` (the integrated variable) from the tide raised on `deformed`.
    /// Leconte et al. (2010) Eq. 2, one of its two symmetric terms:
    /// ```text
    /// (1/a) da/dt = 4 a / (G M_star M_p) K [N(e) omega / n - N_a(e)]
    /// ```
    pub(crate) fn semi_major_axis_13_div_2_derivative(
        &self,
        deformed: &impl ParticleT,
        perturber: &impl ParticleT,
        orbit: &Planet,
    ) -> f64 {
        let e = orbit.eccentricity;
        let relative_rate = 4. * orbit.semi_major_axis
            / (GRAVITATIONAL * perturber.mass() * deformed.mass())
            * self.prefactor(deformed, perturber, orbit)
            * (Self::big_n(e) * deformed.spin() / orbit.mean_motion - Self::big_n_a(e));

        6.5 * orbit.semi_major_axis.powf(6.5) * relative_rate
    }

    /// Derivative of `e^2` (the integrated variable) from the tide raised on `deformed`.
    /// Leconte et al. (2010) Eq. 6, one of its two symmetric terms:
    /// ```text
    /// (1/e) de/dt = 11 a / (G M_star M_p) K [Omega_e(e) omega / n - 18/11 N_e(e)]
    /// ```
    /// so `d(e^2)/dt = 2 e^2 (1/e) de/dt`, which is 0 for a circular orbit.
    pub(crate) fn eccentricity_squared_derivative(
        &self,
        deformed: &impl ParticleT,
        perturber: &impl ParticleT,
        orbit: &Planet,
    ) -> f64 {
        let e = orbit.eccentricity;
        if e == 0. {
            return 0.;
        }
        let relative_rate = 11. * orbit.semi_major_axis
            / (GRAVITATIONAL * perturber.mass() * deformed.mass())
            * self.prefactor(deformed, perturber, orbit)
            * (Self::big_omega_e(e) * deformed.spin() / orbit.mean_motion
                - 18. / 11. * Self::big_n_e(e));

        2. * e * e * relative_rate
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Bolmont et al. 2015, https://doi.org/10.1051/0004-6361/201525909
// Bolmont & Mathis 2016, https://doi.org/10.1007/s10569-016-9690-3
// Hansen 2010, https://doi.org/10.1088/0004-637X/723/1/285
// Hut 1981, https://ui.adsabs.harvard.edu/abs/1981A%26A....99..126H/abstract
// Leconte et al. 2010, https://doi.org/10.1051/0004-6361/201014337
