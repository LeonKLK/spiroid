use crate::constants::MAGNETIC_PERMEABILITY_OF_VACUUM;
use crate::universe::effects::tides::{DynamicalTide, EquilibriumTide, PlanetarySpin, Tides};
use crate::universe::effects::{GeneralRelativityModel, MagneticModel, WindModel};
pub(crate) mod planet;
pub(crate) mod star;

pub use planet::Planet;
pub use star::{Star, StarCsv};

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum ParticleType {
    Planet(Planet),
    Star(Star),
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub struct Particle {
    pub kind: ParticleType,
    #[serde(default)]
    pub tides: Tides,
    #[serde(default)]
    pub(crate) magnetism: MagneticModel,
    #[serde(default)]
    pub(crate) wind: WindModel,
    #[serde(default)]
    pub(crate) general_relativity: GeneralRelativityModel,
}

impl Particle {
    pub(crate) fn is_star(&self) -> bool {
        matches!(self.kind, ParticleType::Star(_))
    }

    pub(crate) fn is_planet(&self) -> bool {
        matches!(self.kind, ParticleType::Planet(_))
    }

    pub(crate) fn initialise(&mut self, time: f64) -> Result<()> {
        match &mut self.kind {
            ParticleType::Star(star) => {
                // Missing inputs deserialise to 0 (serde default): refuse to run the wind with
                // an unset prefactor instead of silently applying no braking.
                if self.wind.wind_torque() && !(star.wind_torque_prefactor > 0.0) {
                    bail!(
                        "wind is enabled but the star input `wind_torque_prefactor` (J) is missing or not positive"
                    );
                }
                if let Some(ctl) = self.tides.constant_time_lag_get() {
                    // `k2_delta_t` is a lag constant in time: for the star the lag follows the
                    // radius through `sigma_bar_star` (Hansen 2010, Bolmont et al. 2015 Eq. 28).
                    if ctl.k2_delta_t != 0.0 {
                        bail!(
                            "`k2_delta_t` is the planetary input: the constant time lag tide of the star takes `sigma_bar_star`"
                        );
                    }
                    // Missing inputs deserialise to 0 (serde default); NaN is refused too.
                    if ctl.sigma_bar_star.partial_cmp(&0.0) != Some(Ordering::Greater) {
                        bail!(
                            "the stellar constant time lag tide needs a positive `sigma_bar_star`"
                        );
                    }
                }
                star.initialise(time)?;
            }
            ParticleType::Planet(planet) => {
                if matches!(self.tides.equilibrium, EquilibriumTide::Zahn(_)) {
                    bail!(
                        "the `Zahn` equilibrium tide is a stellar tide: set it on the central body; the planet takes `ConstantTimeLag`"
                    );
                }
                if matches!(self.tides.dynamical, DynamicalTide::FrequencyAveraged) {
                    bail!(
                        "the `FrequencyAveraged` dynamical tide is a stellar tide: set it on the central body; the planet takes `KaulaTides`"
                    );
                }
                if let Some(ctl) = self.tides.constant_time_lag_get() {
                    // The love number of the kaula planetary tide already holds the whole
                    // response of the planet, and it integrates the planet spin.
                    if self.tides.kaula_enabled() {
                        bail!(
                            "one tidal model on the planet: `ConstantTimeLag` and `KaulaTides` cannot both be enabled on the orbiting body"
                        );
                    }
                    // `sigma_bar_star` is normalised with solar values (Hansen 2010 Eq. 7).
                    if ctl.sigma_bar_star != 0.0 {
                        bail!(
                            "`sigma_bar_star` is the stellar dissipation: the constant time lag tide of the planet takes `k2_delta_t` (s)"
                        );
                    }
                    // Missing inputs deserialise to 0 (serde default); NaN is refused too.
                    if ctl.k2_delta_t.partial_cmp(&0.0) != Some(Ordering::Greater) {
                        bail!(
                            "the planetary constant time lag tide needs a positive `k2_delta_t` (s)"
                        );
                    }
                    // A free spin is evolved with the planet's moment of inertia.
                    if ctl.spin == PlanetarySpin::Free
                        && planet.radius_of_gyration_2.partial_cmp(&0.0) != Some(Ordering::Greater)
                    {
                        bail!(
                            "the planetary constant time lag tide with a `Free` spin needs the planet input `radius_of_gyration_2` (0.254 for Jupiter) to be positive"
                        );
                    }
                }
                planet.initialise();
            }
        }

        Ok(())
    }
}

// Common properties of both Star and Planet.
// Enables making functions generic over impl ParticleT.
pub trait ParticleT {
    fn semi_major_axis(&self) -> f64;
    fn mass(&self) -> f64;
    fn radius(&self) -> f64;
    fn spin(&self) -> f64;
    fn spin_inclination(&self) -> f64;
    fn eccentricity(&self) -> f64;
    fn inclination(&self) -> f64;
    fn luminosity(&self) -> f64;
    fn mean_motion(&self) -> f64;
    fn moment_of_inertia(&self) -> f64;
    fn reduced_mass(&self) -> f64;
}

// https://en.wikipedia.org/wiki/Magnetic_pressure
pub(crate) fn magnetic_pressure(magnetic_field: f64) -> f64 {
    magnetic_field.powi(2) / (2. * MAGNETIC_PERMEABILITY_OF_VACUUM)
}
