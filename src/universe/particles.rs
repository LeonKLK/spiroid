use crate::constants::MAGNETIC_PERMEABILITY_OF_VACUUM;
use crate::universe::effects::tides::{PlanetarySpin, TidalModel};
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
    pub tides: TidalModel,
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
                if self.tides.planetary_constant_time_lag_enabled() {
                    bail!(
                        "`PlanetaryConstantTimeLag` is a planetary tide: set it on the orbiting body; the star takes `StellarConstantTimeLag` or `KaulaTides`"
                    );
                }
                star.initialise(time)?;
            }
            ParticleType::Planet(planet) => {
                if let Some(ctl) = self.tides.planetary_constant_time_lag_get() {
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
