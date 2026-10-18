use anyhow::Result;
use num_complex::Complex;
use sci_file::{DataStore, read_json_from_file};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
mod love_number;
mod polynomials;

pub use love_number::SpectrumFile;
use love_number::{LoveNumber, ParticleComposition, ThermalTideAtmosphereModel};
use polynomials::Polynomials;

use crate::constants::GRAVITATIONAL;
use crate::universe::particles::{ParticleT, Planet};
use derive_more::Add;

// Upper and lower bound for the m, p, q summation.
// Calculated at each time step, based on inclination and eccentricity.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone, Copy)]
pub(crate) struct Mpq {
    m_min: u8,
    m_max: u8,
    p_min: u8,
    p_max: u8,
    q_min: u8,
    q_max: u8,
}

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone, Copy, Add)]
struct Summation {
    #[serde(skip)]
    real_2pq_2mp_dt: f64,
    #[serde(skip)]
    real_2pq_dt_2mp: f64,
    #[serde(skip)]
    imaginary_mfactor: f64,
    #[serde(skip)]
    imaginary_pfactor: f64,
    #[serde(skip)]
    imaginary_qfactor: f64,
    #[serde(skip)]
    imaginary_eccentricity: f64,
    #[serde(skip)]
    imaginary_inclination: f64,
}

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
#[serde(default)]
#[serde(deny_unknown_fields)]
pub struct Kaula {
    pub(crate) particle_type: ParticleComposition,
    pub(crate) atmosphere_model: ThermalTideAtmosphereModel,
    // Inclination and eccentricty polynomials
    #[serde(skip)]
    polynomials: Polynomials,
    #[serde(skip)]
    love_number: LoveNumber,
    // current cache
    #[serde(skip)]
    summation: Summation,

    // cached previous tidal deformed body values
    // Calculations depending on any of the below values are only recomputed when
    // the values have changed across timesteps.
    #[serde(skip)]
    prev_mean_motion: f64,
    #[serde(skip)]
    prev_spin: f64,
    #[serde(skip)]
    prev_inclination: f64,
    #[serde(skip)]
    prev_eccentricity: f64,
}

impl Kaula {
    pub fn interpolation_mode(&self) -> bool {
        !matches!(self.particle_type, ParticleComposition::None)
    }

    pub fn solid_file(&mut self) -> Option<(&PathBuf, &mut DataStore<Complex<f64>>)> {
        match self.particle_type {
            ParticleComposition::Solid {
                ref solid_file,
                ref mut solid_k2,
                ..
            } => Some((solid_file, solid_k2)),
            _ => None,
        }
    }

    pub fn liquid_file(&mut self) -> Option<(&PathBuf, &mut DataStore<Complex<f64>>)> {
        match self.particle_type {
            ParticleComposition::Liquid {
                ref liquid_file,
                ref mut liquid_k2,
                ..
            } => Some((liquid_file, liquid_k2)),
            _ => None,
        }
    }

    pub fn spectrum_file(&mut self) -> Option<(&PathBuf, &mut f64, &mut DataStore<Complex<f64>>)> {
        match self.particle_type {
            ParticleComposition::StarConvectiveEnvelope {
                ref spectrum_file,
                ref mut spin_spec,
                ref mut spectrum_k2,
            } => Some((spectrum_file, spin_spec, spectrum_k2)),
            _ => None,
        }
    }

    /// Loads the love number spectrum of `ParticleComposition::StarConvectiveEnvelope`
    /// from its `SpectrumFile` (no-op for the other compositions).
    /// The file stores dimensionless tidal frequencies (omega / `spin_spec`); they are
    /// converted to rad.s-1 at the reference spin, as expected by `LoveNumber::compute_k2`.
    pub fn load_spectrum_file(&mut self) -> Result<()> {
        if let Some((spectrum_file, spin_spec, spectrum_k2)) = self.spectrum_file() {
            let file: SpectrumFile = read_json_from_file(spectrum_file)?;
            *spin_spec = file.spin_spec;
            *spectrum_k2 = file.spectrum;
            if let DataStore::Interpolate1D(interpolator) = spectrum_k2 {
                interpolator
                    .x_vals_mut()
                    .iter_mut()
                    .for_each(|x| *x *= file.spin_spec);
            }
            spectrum_k2.dimension_check()?;
        }
        Ok(())
    }

    /// Tidal torque exerted on the tidally deformed body by the perturber (kg.m2.s-2).
    /// Boue & Efroimksy (2019) Eq. 123 and Revol et al. (2023) Eq A.3, i.e. the torque
    /// whose division by the moment of inertia gives the spin derivative of the deformed
    /// body (see `planet_spin_derivative` in physics.rs for the planet-role version).
    /// For the stellar tide, pass the star as `tidal_deformed_body` and the planet as
    /// `tidal_perturber` and `orbit`; the result feeds `Star::tidal_torque_convective`.
    pub(crate) fn tidal_torque_on_deformed_body(
        &self,
        tidal_deformed_body: &impl ParticleT,
        tidal_perturber: &impl ParticleT,
        orbit: &Planet,
    ) -> f64 {
        (GRAVITATIONAL * tidal_perturber.mass().powi(2) * tidal_deformed_body.radius().powi(5))
            / orbit.semi_major_axis().powi(6)
            * self.summation_of_longitudinal_modes_spin()
    }

    /// Roles (see `refresh`): `tidal_deformed_body` provides body quantities,
    /// `tidal_perturber` the perturbing mass, `orbit` the orbital elements.
    pub fn initialise_cache(
        &mut self,
        time: f64,
        tidal_perturber: &impl ParticleT,
        tidal_deformed_body: &impl ParticleT,
        orbit: &Planet,
    ) -> Result<()> {
        // Initialise to unreachable values, which forces the caches to be calculated here,
        // reusing the `refresh` function
        // Because comparison with NAN always returns false, see IEEE 754.
        self.prev_eccentricity = f64::NAN;
        self.prev_inclination = f64::NAN;
        self.prev_spin = f64::NAN;
        self.prev_mean_motion = f64::NAN;
        self.refresh(time, tidal_deformed_body, tidal_perturber, orbit)?;
        // No longer NAN here...

        Ok(())
    }

    fn bound_q_by_eccentricity(eccentricity: f64) -> (u8, u8) {
        match () {
            // Select the order of the summation q over the eccentricity function G_lpq
            () if (eccentricity > 0.30) => (0, 15), // q: -7 <= q <= 7
            () if (eccentricity > 0.25) => (1, 14), // q: -6 <= q <= 6
            () if (eccentricity > 0.20) => (2, 13), // q: -5 <= q <= 5
            () if (eccentricity > 0.15) => (3, 12), // q: -4 <= q <= 4
            () if (eccentricity > 1e-8) => (5, 10), // q: -2 <= q <= 2
            () if (eccentricity > 0.0) => (6, 9),   // q: -1 <= q <= 1
            () if (eccentricity == 0.0) => (7, 8),  // q:  0 <= q <= 0
            () => unreachable!("eccentricity cannot be negative."),
        }
    }

    #[allow(clippy::float_cmp)]
    fn eccentricity_changed(&self, orbit: &Planet) -> bool {
        self.prev_eccentricity != orbit.eccentricity()
    }

    #[allow(clippy::float_cmp)]
    fn inclination_changed(&self, orbit: &Planet) -> bool {
        self.prev_inclination != orbit.inclination()
    }

    fn refresh_polynomials(&mut self, orbit: &Planet) {
        if self.eccentricity_changed(orbit) {
            self.polynomials
                .refresh_eccentricity_cache(orbit.eccentricity());
        }
        if self.inclination_changed(orbit) {
            self.polynomials
                .refresh_inclination_cache(orbit.inclination());
        }
    }

    // Save the parameters for comparison during the next time point.
    fn save_parameters(&mut self, tidal_deformed_body: &impl ParticleT, orbit: &Planet) {
        self.prev_eccentricity = orbit.eccentricity();
        self.prev_inclination = orbit.inclination();
        self.prev_mean_motion = orbit.mean_motion();
        self.prev_spin = tidal_deformed_body.spin();
    }

    // Only recalculate if any of the values used in the computation of k2 changed.
    #[allow(clippy::float_cmp)]
    fn love_number_recalculation_needed(
        &self,
        tidal_deformed_body: &impl ParticleT,
        orbit: &Planet,
    ) -> bool {
        self.prev_mean_motion != orbit.mean_motion() || self.prev_spin != tidal_deformed_body.spin()
    }

    // Only recalculate if any of the values used in the computation of the polynomials changed.
    fn summation_recalculation_needed(
        &self,
        tidal_deformed_body: &impl ParticleT,
        orbit: &Planet,
    ) -> bool {
        self.eccentricity_changed(orbit)
            || self.inclination_changed(orbit)
            || self.love_number_recalculation_needed(tidal_deformed_body, orbit)
    }

    fn refresh_summation(
        &mut self,
        time: f64,
        tidal_deformed_body: &impl ParticleT,
        tidal_perturber: &impl ParticleT,
        orbit: &Planet,
        mpq: Mpq,
        summation: &mut Summation,
    ) -> Result<()> {
        // Only recalculate if any of the values used in the computation of k2 changed.
        if self.love_number_recalculation_needed(tidal_deformed_body, orbit) {
            self.love_number.refresh_cache(
                time,
                tidal_deformed_body,
                tidal_perturber,
                orbit,
                &self.particle_type,
                &self.atmosphere_model,
                mpq,
            )?;
        }

        // Only recalculate if inclination or eccentricity changed.
        if self.summation_recalculation_needed(tidal_deformed_body, orbit) {
            summation.imaginary_mfactor = self.sum_over_m_imaginary_mfactor(mpq);
            summation.imaginary_qfactor = self.sum_over_m_imaginary_qfactor(mpq);

            if orbit.eccentricity() != 0.0 {
                summation.imaginary_eccentricity =
                    self.sum_over_m_imaginary_eccentricity(orbit, mpq);
                summation.real_2pq_dt_2mp = self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared_derivative,
                    &self.polynomials.inclination_2mp_squared,
                    mpq,
                );
            }

            if orbit.inclination() != 0.0 && tidal_deformed_body.spin_inclination() != 0.0 {
                summation.imaginary_pfactor = self.sum_over_m_imaginary_pfactor(mpq);
                summation.real_2pq_2mp_dt = self.sum_over_m_real(
                    &self.polynomials.eccentricity_2pq_squared,
                    &self.polynomials.inclination_2mp_squared_derivative,
                    mpq,
                );
            }
        }

        if orbit.inclination() != 0.0 && sin!(orbit.inclination()) != 0.0 {
            summation.imaginary_inclination =
                self.sum_over_m_imaginary_inclination(tidal_deformed_body, orbit, mpq);
        }
        Ok(())
    }

    // All the calculations using the polynomials and love number are performed here
    // and stored in the `sum_over_xxx` caches.
    // The caches are used when the derivitaves are calculated for each keplerian element.
    //
    // Three roles:
    // - `tidal_deformed_body`: the body raising the tidal bulge (planet for the planetary
    //   tide, star for the stellar tide). Provides spin, radius, mass, moment of inertia
    //   and spin inclination.
    // - `tidal_perturber`: the body raising the tide (mass, luminosity).
    // - `orbit`: the body carrying the orbital elements (mean motion, semi-major axis,
    //   eccentricity, inclination, reduced mass). Always the planet, since the orbit is
    //   stored on the planet regardless of which body is deformed.
    // For the planetary tide `tidal_deformed_body` and `orbit` are the same planet.
    pub(crate) fn refresh(
        &mut self,
        time: f64,
        tidal_deformed_body: &impl ParticleT,
        tidal_perturber: &impl ParticleT,
        orbit: &Planet,
    ) -> Result<()> {
        self.refresh_polynomials(orbit);
        let (q_min, q_max) = Self::bound_q_by_eccentricity(orbit.eccentricity());

        if orbit.inclination() <= 1e-8 {
            // If inclination is close to zero, only compute
            // m = 0, p = 1 and m = 2, p = 0
            let mpq_01q = Mpq {
                m_min: 0,
                m_max: 1,
                p_min: 1,
                p_max: 2,
                q_min,
                q_max,
            };

            let mpq_20q = Mpq {
                m_min: 2,
                m_max: 3,
                p_min: 0,
                p_max: 1,
                q_min,
                q_max,
            };

            let mut summation_01q = self.summation;
            self.refresh_summation(
                time,
                tidal_deformed_body,
                tidal_perturber,
                orbit,
                mpq_01q,
                &mut summation_01q,
            )?;
            // Start with zero values, any that require updating will be updated
            let mut summation_20q = Summation::default();
            self.refresh_summation(
                time,
                tidal_deformed_body,
                tidal_perturber,
                orbit,
                mpq_20q,
                &mut summation_20q,
            )?;
            // Add the two summations. Non-updated values will contain originals in 01q and 0 in 20q.
            self.summation = summation_01q + summation_20q;
        } else {
            let mpq = Mpq {
                m_min: 0,
                m_max: 3,
                p_min: 0,
                p_max: 3,
                q_min,
                q_max,
            };

            let mut summation = self.summation;
            self.refresh_summation(
                time,
                tidal_deformed_body,
                tidal_perturber,
                orbit,
                mpq,
                &mut summation,
            )?;
            self.summation = summation;
        }
        //        panic!();
        self.save_parameters(tidal_deformed_body, orbit);
        Ok(())
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the semi-major-axis derivative.
    // Boue & Efroimksy (2019) Eq 116 and Revol et al. (2023) Eq A.1.
    pub(crate) fn summation_of_longitudinal_modes_semi_major_axis(&self) -> f64 {
        self.summation.imaginary_qfactor
    }

    // Summation over longitudinal modes m for the computation of the spin derivative.
    // Boue & Efroimksy (2019) Eq 123 and Revol et al. (2023) Eq A.3
    pub(crate) fn summation_of_longitudinal_modes_spin(&self) -> f64 {
        self.summation.imaginary_mfactor
    }

    pub(crate) fn summation_of_longitudinal_modes_eccentricity(&self) -> f64 {
        self.summation.imaginary_eccentricity
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the eccentricity derivative.
    // Boue & Efroimksy (2019) Eq 117 and Revol et al. (2023) Eq A.3
    fn sum_over_m_imaginary_eccentricity(&self, orbit: &Planet, mpq: Mpq) -> f64 {
        let semi_minor_axis_ratio = sqrt!(1. - orbit.eccentricity().powi(2));

        self.polynomials
            .inclination_2mp_squared
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
            .map(|(m, m_val)| {
                self.polynomials
                    .eccentricity_2pq_squared
                    .iter()
                    .zip(m_val)
                    .enumerate()
                    .take(mpq.p_max.into())
                    .skip(mpq.p_min.into())
                    .map(|(p, (p_val, m_val_p))| {
                        let p_factor = (2 - 2 * (p as isize)) as f64;
                        p_val
                            .iter()
                            .enumerate()
                            .take(mpq.q_max.into())
                            .skip(mpq.q_min.into())
                            .map(|(q, q_val)| {
                                let q_factor = p_factor + (q as f64 - 7.);
                                let term = q_factor * semi_minor_axis_ratio - p_factor;
                                let imk2 = self.love_number.k2(m, p, q).im;
                                imk2 * q_val * term
                            })
                            .sum::<f64>()
                            * m_val_p
                    })
                    .sum::<f64>()
                    * factorial_kronecker(m)
            })
            .sum::<f64>()
    }

    pub(crate) fn summation_of_longitudinal_modes_inclination(&self) -> f64 {
        self.summation.imaginary_inclination
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the inclination derivative.
    // by Boue & Efroimksy (2019) Eq 118 and Revol et al. (2023) Eq A.7
    fn sum_over_m_imaginary_inclination(
        &self,
        tidal_deformed_body: &impl ParticleT,
        orbit: &Planet,
        mpq: Mpq,
    ) -> f64 {
        let semi_minor_axis_ratio = sqrt!(1. - orbit.eccentricity().powi(2));

        let term1 =
            (orbit.reduced_mass() * orbit.mean_motion().powi(2) * orbit.semi_major_axis().powi(2))
                / (tidal_deformed_body.moment_of_inertia() * tidal_deformed_body.spin());
        let term3 = orbit.mean_motion() / semi_minor_axis_ratio;
        let cos_inc = cos!(orbit.inclination());

        self.polynomials
            .inclination_2mp_squared
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
            .map(|(m, m_val)| {
                self.polynomials
                    .eccentricity_2pq_squared
                    .iter()
                    .zip(m_val)
                    .enumerate()
                    .take(mpq.p_max.into())
                    .skip(mpq.p_min.into())
                    .map(|(p, (p_val, m_val_p))| {
                        let p_factor = (2 - 2 * (p as isize)) as f64;
                        p_val
                            .iter()
                            .enumerate()
                            .take(mpq.q_max.into())
                            .skip(mpq.q_min.into())
                            .map(|(q, q_val)| self.love_number.k2(m, p, q).im * q_val)
                            .sum::<f64>()
                            * m_val_p
                            * (term1 * (m as f64 * cos_inc - p_factor)
                                - ((p_factor * cos_inc - m as f64) * term3))
                    })
                    .sum::<f64>()
                    * factorial_kronecker(m)
            })
            .sum::<f64>()
    }

    // TODO (stellar tide, non-coplanar case):
    // Three of the summations below (longitude of ascending node, spin axis inclination,
    // pericentre inclination) take a `&Planet` and read from it BOTH the
    // deformed-body quantities (`moment_of_inertia`, `spin`, `tan_spin_inc`) AND the
    // orbital quantities (`tan_inc`, `sin_inc`, `sin_lan`, `cos_lan`, `mean_motion`,
    // `semi_major_axis`, `semi_minor_axis_ratio`, `reduced_mass`). This is only correct
    // when the planet is the tidally deformed body. When the star is the deformed body,
    // the spin-axis terms refer to the star, which has no `spin_inclination` /
    // `longitude_ascending_node` state (its spin axis is the reference frame, see
    // `Star::spin_inclination`), so `1 / (I * spin * tan_spin_inc)` is not defined as-is.
    // These functions must be split into (tidal_deformed_body, orbit) roles, like
    // `sum_over_m_imaginary_inclination`, once the non-coplanar stellar tide is needed.
    //
    // For now the coplanar case is sufficient: with `inclination == 0` and
    // `spin_inclination == 0` none of these three functions is ever called (see the gates
    // in `physics.rs`). The summations used by the stellar tide are
    // `summation_of_longitudinal_modes_{semi_major_axis, spin, eccentricity}`, which are
    // argument-free, and `summation_of_longitudinal_modes_pericentre_eccentricity`, which
    // only reads orbital quantities; all four are role-agnostic.
    fn summation_of_longitudinal_modes_triple_common(
        &self,
        term1: f64,
        term2: f64,
        term3: f64,
    ) -> f64 {
        (self.summation.real_2pq_2mp_dt * term1 * 0.5)
            + (self.summation.imaginary_mfactor * term2)
            + (self.summation.imaginary_pfactor * term3)
    }

    // Summation over longitudinal modes m for the computation of the longitude of ascending node derivative.
    // Boue & Efroimksy (2019) Eq 121 and Revol et al. (2023) Eq A.9
    pub(crate) fn summation_of_longitudinal_modes_longitude_ascending_node(
        &self,
        tidal_deformed_body: &Planet,
    ) -> f64 {
        let term1 = (1.
            / (tidal_deformed_body.moment_of_inertia
                * tidal_deformed_body.spin
                * tidal_deformed_body.tan_inc))
            - (tidal_deformed_body.cos_lan
                / (tidal_deformed_body.moment_of_inertia
                    * tidal_deformed_body.spin
                    * tidal_deformed_body.tan_spin_inc))
            + (1.
                / (tidal_deformed_body.reduced_mass
                    * tidal_deformed_body.mean_motion
                    * tidal_deformed_body.semi_major_axis.powi(2)
                    * tidal_deformed_body.semi_minor_axis_ratio
                    * tidal_deformed_body.sin_inc));

        let term2 = -(tidal_deformed_body.sin_lan * cotan!(tidal_deformed_body.inclination))
            / (tidal_deformed_body.moment_of_inertia
                * tidal_deformed_body.spin
                * tidal_deformed_body.tan_spin_inc);

        let term3 = tidal_deformed_body.sin_lan
            / (tidal_deformed_body.moment_of_inertia
                * tidal_deformed_body.spin
                * tidal_deformed_body.tan_spin_inc
                * tidal_deformed_body.sin_inc);

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the spin axis inclination derivative.
    // Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.12
    pub(crate) fn summation_of_longitudinal_modes_spin_axis_inclination(
        &self,
        tidal_deformed_body: &Planet,
    ) -> f64 {
        let term1 = -tidal_deformed_body.sin_lan;
        let term2 = tidal_deformed_body.cos_lan / tidal_deformed_body.tan_inc;
        let term3 = -(tidal_deformed_body.cos_lan / tidal_deformed_body.sin_inc);

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the eccentricity dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    // Only orbital quantities are read, so this is valid for both the planetary tide
    // (planet deformed) and the stellar tide (star deformed); `orbit` is always the planet.
    pub(crate) fn summation_of_longitudinal_modes_pericentre_eccentricity(
        &self,
        orbit: &Planet,
    ) -> f64 {
        let term2 = orbit.semi_minor_axis_ratio
            / (orbit.mean_motion
                * orbit.semi_major_axis.powi(2)
                * orbit.eccentricity
                * orbit.reduced_mass);
        self.summation.real_2pq_dt_2mp * 0.5 * term2
    }

    // Summation over longitudinal modes m for the computation of the inclination dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_inclination(
        &self,
        tidal_deformed_body: &Planet,
    ) -> f64 {
        let term1 = -((1.
            / (tidal_deformed_body.moment_of_inertia
                * tidal_deformed_body.spin
                * tidal_deformed_body.sin_inc))
            + (1.
                / (tidal_deformed_body.mean_motion
                    * tidal_deformed_body.semi_major_axis.powi(2)
                    * tidal_deformed_body.semi_minor_axis_ratio
                    * tidal_deformed_body.tan_inc
                    * tidal_deformed_body.reduced_mass)));
        self.summation.real_2pq_2mp_dt * 0.5 * term1
    }

    // Iteration over the provided 2D arrays (outer 3x15 and inner 3x3), summing the contents of:
    // (love_number(m, p, q) * inner[p][q]) * outer[m][p] * factorial_kronecker(m)
    fn sum_over_m_real(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3], mpq: Mpq) -> f64 {
        outer
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
            .map(|(m, m_val)| {
                inner
                    .iter()
                    .zip(m_val)
                    .enumerate()
                    .take(mpq.p_max.into())
                    .skip(mpq.p_min.into())
                    .map(|(p, (p_val, m_val_p))| {
                        p_val
                            .iter()
                            .enumerate()
                            .take(mpq.q_max.into())
                            .skip(mpq.q_min.into())
                            .map(|(q, q_val)| self.love_number.k2(m, p, q).re * q_val)
                            .sum::<f64>()
                            * m_val_p
                    })
                    .sum::<f64>()
                    * factorial_kronecker(m)
            })
            .sum::<f64>()
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    fn sum_over_m_imaginary_pfactor(&self, mpq: Mpq) -> f64 {
        self.polynomials
            .inclination_2mp_squared
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
            .map(|(m, m_val)| {
                self.polynomials
                    .eccentricity_2pq_squared
                    .iter()
                    .zip(m_val)
                    .enumerate()
                    .take(mpq.p_max.into())
                    .skip(mpq.p_min.into())
                    .map(|(p, (p_val, m_val_p))| {
                        let p_factor = (2 - 2 * (p as isize)) as f64;
                        p_val
                            .iter()
                            .enumerate()
                            .take(mpq.q_max.into())
                            .skip(mpq.q_min.into())
                            .map(|(q, q_val)| self.love_number.k2(m, p, q).im * q_val * p_factor)
                            .sum::<f64>()
                            * m_val_p
                    })
                    .sum::<f64>()
                    * factorial_kronecker(m)
            })
            .sum::<f64>()
    }

    // Precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    fn sum_over_m_imaginary_mfactor(&self, mpq: Mpq) -> f64 {
        // Skip over the case of m = 0, since it would be 0.
        let m_min = max!(1, mpq.m_min);

        self.polynomials
            .inclination_2mp_squared
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(m_min.into())
            .map(|(m, m_val)| {
                self.polynomials
                    .eccentricity_2pq_squared
                    .iter()
                    .zip(m_val)
                    .enumerate()
                    .take(mpq.p_max.into())
                    .skip(mpq.p_min.into())
                    .map(|(p, (p_val, m_val_p))| {
                        p_val
                            .iter()
                            .enumerate()
                            .take(mpq.q_max.into())
                            .skip(mpq.q_min.into())
                            .map(|(q, q_val)| self.love_number.k2(m, p, q).im * q_val)
                            .sum::<f64>()
                            * m_val_p
                    })
                    .sum::<f64>()
                    * factorial_kronecker(m)
                    * (m as f64)
            })
            .sum::<f64>()
    }

    // Precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    fn sum_over_m_imaginary_qfactor(&self, mpq: Mpq) -> f64 {
        self.polynomials
            .inclination_2mp_squared
            .iter()
            .enumerate()
            .take(mpq.m_max.into())
            .skip(mpq.m_min.into())
            .map(|(m, m_val)| {
                self.polynomials
                    .eccentricity_2pq_squared
                    .iter()
                    .zip(m_val)
                    .enumerate()
                    .take(mpq.p_max.into())
                    .skip(mpq.p_min.into())
                    .map(|(p, (p_val, m_val_p))| {
                        let p_factor = (2 - 2 * (p as isize)) as f64;
                        p_val
                            .iter()
                            .enumerate()
                            .take(mpq.q_max.into())
                            .skip(mpq.q_min.into())
                            .map(|(q, q_val)| {
                                let q_factor = p_factor + (q as f64 - 7.);
                                self.love_number.k2(m, p, q).im * q_val * q_factor
                            })
                            .sum::<f64>()
                            * m_val_p
                    })
                    .sum::<f64>()
                    * factorial_kronecker(m)
            })
            .sum::<f64>()
    }
}

// This is a precomputed simplification of the original calculation:
//      (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0))
// where factorial(z) is defined recursively as:
//      if z == 0 -> 1
//      else -> z * factorial(z - 1)
// where kronecker_delta(x, y) is defined as:
//      x == y -> 1
//      x != y -> 0
fn factorial_kronecker(m: usize) -> f64 {
    match m {
        0 => 1.,
        1 => 1. / 3.,
        2 => 1. / 12.,
        _ => unreachable!(),
    }
}

#[cfg(test)]
pub mod tests;
