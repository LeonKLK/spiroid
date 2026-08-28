pub(crate) mod kaula;
pub use kaula::{Kaula, SpectrumFile};

pub(crate) mod constant_time_lag;
pub use constant_time_lag::ConstantTimeLag;

use crate::universe::particles::{Planet, Star};

use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
#[allow(clippy::large_enum_variant)]
pub enum TidalModel {
    #[default]
    Disabled,
    // Equilibrium tide dissipation given as the dimensionless sigma_bar_star from Bolmont & Mathis (2016), Eq. 8
    ConstantTimeLag(ConstantTimeLag),
    KaulaTides(Kaula),
}

impl TidalModel {
    pub(crate) fn stellar_tidal_torque(&self, star: &Star, planet: &Planet) -> f64 {
        match self {
            TidalModel::Disabled => 0.0,
            TidalModel::ConstantTimeLag(constant_time_lag) => {
                // requires tidal_frequency
                constant_time_lag.stellar_tidal_torque(star, planet)
            }
            // Stellar tide with the kaula model: the star is the deformed body.
            // Requires `refresh_kaula_star` to have been called for the current state.
            TidalModel::KaulaTides(kaula) => {
                kaula.tidal_torque_on_deformed_body(star, planet, planet)
            }
        }
    }

    /// Refreshes the kaula tides data (love number, eccentricity and inclination polynomials)
    /// for the planetary tide: the planet is the tidally deformed body, the star the
    /// perturber, and the planet also carries the orbit.
    pub(crate) fn refresh_kaula(&mut self, time: f64, star: &Star, planet: &Planet) -> Result<()> {
        if let &mut TidalModel::KaulaTides(ref mut kaula) = self {
            kaula.refresh(time, planet, star, planet)?;
        }

        Ok(())
    }

    /// Refreshes the kaula tides data for the stellar tide: the star is the tidally
    /// deformed body, the planet the perturber, and the planet carries the orbit.
    pub(crate) fn refresh_kaula_star(
        &mut self,
        time: f64,
        star: &Star,
        planet: &Planet,
    ) -> Result<()> {
        if let &mut TidalModel::KaulaTides(ref mut kaula) = self {
            kaula.refresh(time, star, planet, planet)?;
        }

        Ok(())
    }

    /// Returns `true` if the `TidalModel` is `KaulaTides`.
    pub(crate) fn kaula_enabled(&self) -> bool {
        matches!(&self, TidalModel::KaulaTides(_))
    }

    /// Returns a reference to the `Kaula` struct if the `TidalModel` is `KaulaTides`.
    pub fn kaula_get(&self) -> Option<&Kaula> {
        match self {
            TidalModel::KaulaTides(kaula) => Some(kaula),
            _ => None,
        }
    }

    /// Returns a mutable reference to the `Kaula` struct if the `TidalModel` is `KaulaTides`.
    pub fn kaula_get_mut(&mut self) -> Option<&mut Kaula> {
        match self {
            &mut TidalModel::KaulaTides(ref mut kaula) => Some(kaula),
            _ => None,
        }
    }
}

// References:
// Bolmont & Mathis 2016, https://doi.org/10.1007/s10569-016-9690-3
