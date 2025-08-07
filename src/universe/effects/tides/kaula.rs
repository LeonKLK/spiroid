use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

mod love_number;
mod polynomials;

use love_number::{LoveNumber, ParticleComposition};
use polynomials::Polynomials;

use crate::universe::particles::{ParticleT, Planet};
use crate::utils::{factorial, kronecker_delta};

use wtf_tides::Layer;

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
#[serde(default)]
#[serde(deny_unknown_fields)]
pub struct Kaula {
    pub(crate) particle_type: ParticleComposition,
    // Inclination and eccentricty polynomials
    #[serde(skip)]
    polynomials: Polynomials,
    #[serde(skip)]
    love_number: LoveNumber,

    //cache
    #[serde(skip)]
    sum_over_m_real_2pq_2mp_dt: f64,
    #[serde(skip)]
    sum_over_m_imaginary_mfactor: f64,
    #[serde(skip)]
    sum_over_m_imaginary_pfactor: f64,
}

impl Kaula {
    pub fn solid_file(&self) -> Option<&PathBuf> {
        match self.particle_type {
            ParticleComposition::Solid { ref solid_file }
            | ParticleComposition::Atmosphere { ref solid_file, .. }
            | ParticleComposition::SolidAtmosphere { ref solid_file, .. }
            | ParticleComposition::SolidOcean { ref solid_file, .. }
            | ParticleComposition::SolidAtmosphereOcean { ref solid_file, .. } => Some(solid_file),
            _ => None,
        }
    }

    pub fn ocean_file(&self) -> Option<&PathBuf> {
        match self.particle_type {
            ParticleComposition::SolidOcean { ref ocean_file, .. }
            | ParticleComposition::SolidAtmosphereOcean { ref ocean_file, .. } => Some(ocean_file),
            _ => None,
        }
    }

    pub fn interpolate_dir(&self) -> Option<&PathBuf> {
        match self.particle_type {
            ParticleComposition::Interpolate {
                ref interpolate_dir,
                ..
            }
            | ParticleComposition::InterpolateAtmosphere {
                ref interpolate_dir,
                ..
            } => Some(interpolate_dir),
            _ => None,
        }
    }

    pub fn initialise_love_number_solid(&mut self, love_solid: &[Vec<f64>]) {
        self.love_number
            .imaginary_solid
            .init(&love_solid[0], &love_solid[1]);
        self.love_number
            .real_solid
            .init(&love_solid[0], &love_solid[2]);
    }

    pub fn initialise_love_number_layers(&mut self, planet_layers: &[Layer]) {
        self.love_number.init_layers(planet_layers);
    }

    pub fn initialise_love_number_ocean(&mut self, love_ocean: &[Vec<f64>]) {
        self.love_number
            .imaginary_oceanic
            .init(&love_ocean[0], &love_ocean[1]);
    }

    pub(crate) fn refresh(
        &mut self,
        time: f64,
        planet: &impl ParticleT,
        star: &impl ParticleT,
    ) -> Result<()> {
        self.polynomials
            .refresh_cache(planet.eccentricity(), planet.inclination());
        self.love_number
            .refresh_cache(time, planet, star, &self.particle_type)?;

        if planet.inclination() != 0.0 {
            self.sum_over_m_real_2pq_2mp_dt = self.sum_over_m_real(
                &self.polynomials.eccentricity_2pq_squared,
                &self.polynomials.inclination_2mp_squared_derivative,
            );
            self.sum_over_m_imaginary_mfactor = self.sum_over_m_imaginary_mfactor(
                &self.polynomials.eccentricity_2pq_squared,
                &self.polynomials.inclination_2mp_squared,
            );
            self.sum_over_m_imaginary_pfactor = self.sum_over_m_imaginary_pfactor(
                &self.polynomials.eccentricity_2pq_squared,
                &self.polynomials.inclination_2mp_squared,
            );
        }

        Ok(())
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the semi-major-axis derivative.
    // Boue & Efroimksy (2019) Eq 116 and Revol et al. (2023) Eq A.1.
    pub(crate) fn summation_of_longitudinal_modes_semi_major_axis(&self) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut sum_over_m = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut sum_over_p = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut sum_over_q = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    let q_factor = p_factor + (q as f64 - 7.);
                    let imk2 = self.love_number.imaginary(m, p, q);
                    sum_over_q += q_factor * q_val * imk2;
                }
                sum_over_p += m_val[p] * sum_over_q;
            }
            sum_over_m +=
                sum_over_p * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        sum_over_m
    }

    // Summation over longitudinal modes m for the computation of the spin derivative.
    // Boue & Efroimksy (2019) Eq 123 and Revol et al. (2023) Eq A.3
    pub(crate) fn summation_of_longitudinal_modes_spin(&self) -> f64 {
        self.sum_over_m_imaginary_mfactor(
            &self.polynomials.eccentricity_2pq_squared,
            &self.polynomials.inclination_2mp_squared,
        )
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the eccentricity derivative.
    // Boue & Efroimksy (2019) Eq 117 and Revol et al. (2023) Eq A.3
    pub(crate) fn summation_of_longitudinal_modes_eccentricity(
        &self,
        semi_minor_axis_ratio: f64,
    ) -> f64 {
        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut sum_over_m = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut sum_over_p = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut sum_over_q = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    let q_factor = p_factor + (q as f64 - 7.);
                    let term = q_factor * semi_minor_axis_ratio - p_factor;
                    let imk2 = self.love_number.imaginary(m, p, q);
                    sum_over_q += imk2 * q_val * term;
                }
                sum_over_p += sum_over_q * m_val[p];
            }
            sum_over_m +=
                sum_over_p * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        sum_over_m
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    // Summation over longitudinal modes m for the computation of the inclination derivative.
    // by Boue & Efroimksy (2019) Eq 118 and Revol et al. (2023) Eq A.7
    pub(crate) fn summation_of_longitudinal_modes_inclination(&self, planet: &Planet) -> f64 {
        let term1 =
            (planet.reduced_mass * planet.mean_motion.powi(2) * planet.semi_major_axis.powi(2))
                / (planet.moment_of_inertia * planet.spin);
        let term3 = planet.mean_motion / planet.semi_minor_axis_ratio;

        let outer = &self.polynomials.inclination_2mp_squared;
        let inner = &self.polynomials.eccentricity_2pq_squared;

        let mut sum_over_m = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut sum_over_p = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut sum_over_q = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    let imk2 = self.love_number.imaginary(m, p, q);
                    sum_over_q += imk2 * q_val;
                }
                sum_over_p += sum_over_q
                    * m_val[p]
                    * (term1 * (m as f64 * planet.cos_inc - p_factor)
                        - ((p_factor * planet.cos_inc - m as f64) * term3));
            }
            sum_over_m +=
                sum_over_p * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        sum_over_m
    }

    fn summation_of_longitudinal_modes_triple_common(
        &self,
        term1: f64,
        term2: f64,
        term3: f64,
    ) -> f64 {
        (self.sum_over_m_real_2pq_2mp_dt * term1 * 0.5)
            + (self.sum_over_m_imaginary_mfactor * term2)
            + (self.sum_over_m_imaginary_pfactor * term3)
    }

    // Summation over longitudinal modes m for the computation of the longitude of ascending node derivative.
    // Boue & Efroimksy (2019) Eq 121 and Revol et al. (2023) Eq A.9
    pub(crate) fn summation_of_longitudinal_modes_longitude_ascending_node(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term1 = (1. / (planet.moment_of_inertia * planet.spin * planet.tan_inc))
            - (planet.cos_lan / (planet.moment_of_inertia * planet.spin * planet.tan_spin_inc))
            + (1.
                / (planet.reduced_mass
                    * planet.mean_motion
                    * planet.semi_major_axis.powi(2)
                    * planet.semi_minor_axis_ratio
                    * planet.sin_inc));

        let term2 = -(planet.sin_lan * cotan!(planet.inclination))
            / (planet.moment_of_inertia * planet.spin * planet.tan_spin_inc);

        let term3 = planet.sin_lan
            / (planet.moment_of_inertia * planet.spin * planet.tan_spin_inc * planet.sin_inc);

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the spin axis inclination derivative.
    // Boue & Efroimksy (2019) Eq 122 and Revol et al. (2023) Eq A.12
    pub(crate) fn summation_of_longitudinal_modes_spin_axis_inclination(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term1 = -planet.sin_lan;
        let term2 = planet.cos_lan / planet.tan_inc;
        let term3 = -(planet.cos_lan / planet.sin_inc);

        self.summation_of_longitudinal_modes_triple_common(term1, term2, term3)
    }

    // Summation over longitudinal modes m for the computation of the eccentricity dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_eccentricity(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term2 = planet.semi_minor_axis_ratio
            / (planet.mean_motion
                * planet.semi_major_axis.powi(2)
                * planet.eccentricity
                * planet.reduced_mass);
        self.sum_over_m_real(
            &self.polynomials.eccentricity_2pq_squared_derivative,
            &self.polynomials.inclination_2mp_squared,
        ) * 0.5
            * term2
    }

    // Summation over longitudinal modes m for the computation of the inclination dependent longitude of pericentre derivative.
    // Boue & Efroimksy (2019) Eq 120 and Revol et al. (2023) Eq A.11
    pub(crate) fn summation_of_longitudinal_modes_pericentre_inclination(
        &self,
        planet: &Planet,
    ) -> f64 {
        let term1 = -((1. / (planet.moment_of_inertia * planet.spin * planet.sin_inc))
            + (1.
                / (planet.mean_motion
                    * planet.semi_major_axis.powi(2)
                    * planet.semi_minor_axis_ratio
                    * planet.tan_inc
                    * planet.reduced_mass)));
        self.sum_over_m_real_2pq_2mp_dt * 0.5 * term1
    }

    // Iteration over the provided 2D arrays (outer 3x15 and inner 3x3), summing the contents of:
    // (love_number(m, p, q) * inner[p][q]) * outer[m]p] * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0))
    fn sum_over_m_real(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3]) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut p_sum = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut q_sum = 0.0;
                for (q, q_val) in p_val.iter().enumerate() {
                    q_sum += q_val * self.love_number.real(m, p, q);
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    // Wrapping and precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    #[allow(clippy::cast_possible_wrap)]
    fn sum_over_m_imaginary_pfactor(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3]) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut p_sum = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut q_sum = 0.0;
                let p_factor = (2 - 2 * (p as isize)) as f64;
                for (q, q_val) in p_val.iter().enumerate() {
                    q_sum += q_val * self.love_number.imaginary(m, p, q) * p_factor;
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum * (factorial(2 - m) / factorial(2 + m)) * (2. - kronecker_delta(m, 0));
        }
        m_sum
    }

    // Precision loss not applicable since the values are in [0..15).
    #[allow(clippy::cast_precision_loss)]
    fn sum_over_m_imaginary_mfactor(&self, inner: &[[f64; 15]; 3], outer: &[[f64; 3]; 3]) -> f64 {
        let mut m_sum = 0.0;
        for (m, m_val) in outer.iter().enumerate() {
            let mut p_sum = 0.0;
            for (p, p_val) in inner.iter().enumerate() {
                let mut q_sum = 0.0;
                for (q, q_val) in p_val.iter().enumerate() {
                    q_sum += q_val * self.love_number.imaginary(m, p, q);
                }
                p_sum += q_sum * m_val[p];
            }
            m_sum += p_sum
                * (factorial(2 - m) / factorial(2 + m))
                * (2. - kronecker_delta(m, 0))
                * (m as f64);
        }
        m_sum
    }
}

#[cfg(test)]
pub mod tests;
