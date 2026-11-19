pub(crate) mod kaula;
pub use kaula::{Kaula, SpectrumFile};

pub(crate) mod constant_time_lag;
pub use constant_time_lag::{ConstantTimeLag, PlanetarySpin};

pub(crate) mod quality_factor;
pub use quality_factor::Zahn;

use crate::universe::particles::{Planet, Star};

use anyhow::Result;
use serde::de::{self, MapAccess, Visitor};
use serde::{Deserialize, Deserializer, Serialize};

/// Tides raised on one body, split by mechanism. Each body carries its own `Tides`: the
/// tide raised on the star is `central_body.tides`, the tide raised on the planet is
/// `orbiting_body.tides`. The two entries are independent and their effects add up.
///
/// In a configuration file:
/// ```json
/// "tides": {
///     "equilibrium": {"ConstantTimeLag": {"sigma_bar_star": 1e-06}},
///     "dynamical": "FrequencyAveraged"
/// }
/// ```
/// A missing entry is disabled, and `"tides": "Disabled"` disables both.
#[derive(Serialize, PartialEq, Debug, Default, Clone)]
pub struct Tides {
    pub equilibrium: EquilibriumTide,
    pub dynamical: DynamicalTide,
}

/// Equilibrium tide: the hydrostatic adjustment of the body to the perturber.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone, Copy)]
pub enum EquilibriumTide {
    #[default]
    Disabled,
    // Constant time lag tide of Hut (1981) / Leconte et al. (2010): any eccentricity,
    // coplanar. The same model for the star and the planet.
    ConstantTimeLag(ConstantTimeLag),
    // Zahn prescription as parameterised in Mustill & Villaver (2012), expressed as a tidal
    // quality factor (Benbakoura et al. 2019): circular orbit. Star only.
    Zahn(Zahn),
}

/// Dynamical tide: the wave-like response of the body.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
#[allow(clippy::large_enum_variant)]
pub enum DynamicalTide {
    #[default]
    Disabled,
    // Frequency-averaged dissipation of inertial waves in the convective envelope (Ogilvie
    // 2013, Mathis 2015), expressed as a tidal quality factor (Bolmont & Mathis 2016,
    // Benbakoura et al. 2019): circular orbit. Star only.
    FrequencyAveraged,
    // Frequency-dependent love number summed over the tidal harmonics (Kaula 1964).
    KaulaTides(Kaula),
}

impl Tides {
    /// Torque exerted on the star by the stellar tides expressed as a tidal quality factor
    /// (`Zahn`, `FrequencyAveraged`). Their dissipations add up (Benbakoura et al. 2019
    /// Eq. 20). On a circular orbit the angular momentum taken from the stellar spin goes
    /// to the orbit, so physics.rs also derives the semi-major axis from this torque.
    pub(crate) fn stellar_quality_factor_torque(&self, star: &Star, planet: &Planet) -> f64 {
        // When a tide is disabled its tidal quality factor would be 1 / 0.
        // We set it to infinity such that 1 / infinity == 0 == disabled.
        let equilibrium_tidal_quality = match &self.equilibrium {
            EquilibriumTide::Zahn(zahn) => zahn.tidal_quality(star, planet),
            EquilibriumTide::Disabled | EquilibriumTide::ConstantTimeLag(_) => f64::INFINITY,
        };
        let dynamical_tidal_quality = match &self.dynamical {
            DynamicalTide::FrequencyAveraged => {
                quality_factor::frequency_averaged::tidal_quality(star)
            }
            DynamicalTide::Disabled | DynamicalTide::KaulaTides(_) => f64::INFINITY,
        };
        if equilibrium_tidal_quality.is_infinite() && dynamical_tidal_quality.is_infinite() {
            return 0.0;
        }
        let total_dissipation = 1. / equilibrium_tidal_quality + 1. / dynamical_tidal_quality;

        // requires tidal_frequency
        quality_factor::stellar_tidal_torque(1. / total_dissipation, star, planet)
    }

    /// Total tidal torque exerted on the star (kg.m2.s-2), applied to its convective zone:
    /// the sum of the enabled stellar tides. Recorded in the output as
    /// `star.tidal_torque_convective`; the integration does not read it back, `physics::force`
    /// collects the same torques in the same order (`stellar_tidal_torque`).
    pub(crate) fn stellar_tidal_torque(&self, star: &Star, planet: &Planet) -> f64 {
        let mut torque = self.stellar_quality_factor_torque(star, planet);

        // Constant time lag equilibrium tide: star deformed, planet perturber.
        if let EquilibriumTide::ConstantTimeLag(ctl) = &self.equilibrium {
            torque += ctl.tidal_torque(star, planet, planet);
        }
        // Stellar tide with the kaula model: the star is the deformed body.
        // Requires `refresh_kaula_star` to have been called for the current state.
        if let DynamicalTide::KaulaTides(kaula) = &self.dynamical {
            torque += kaula.tidal_torque_on_deformed_body(star, planet, planet);
        }

        torque
    }

    /// Refreshes the kaula tides data (love number, eccentricity and inclination polynomials)
    /// for the planetary tide: the planet is the tidally deformed body, the star the
    /// perturber, and the planet also carries the orbit.
    pub(crate) fn refresh_kaula(&mut self, time: f64, star: &Star, planet: &Planet) -> Result<()> {
        if let &mut DynamicalTide::KaulaTides(ref mut kaula) = &mut self.dynamical {
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
        if let &mut DynamicalTide::KaulaTides(ref mut kaula) = &mut self.dynamical {
            kaula.refresh(time, star, planet, planet)?;
        }

        Ok(())
    }

    /// Returns a reference to the `ConstantTimeLag` struct if that is the equilibrium tide.
    pub fn constant_time_lag_get(&self) -> Option<&ConstantTimeLag> {
        match &self.equilibrium {
            EquilibriumTide::ConstantTimeLag(ctl) => Some(ctl),
            _ => None,
        }
    }

    /// Returns `true` if the equilibrium tide is `ConstantTimeLag`.
    pub(crate) fn constant_time_lag_enabled(&self) -> bool {
        self.constant_time_lag_get().is_some()
    }

    /// Returns `true` if the dynamical tide is `KaulaTides`.
    pub(crate) fn kaula_enabled(&self) -> bool {
        matches!(&self.dynamical, DynamicalTide::KaulaTides(_))
    }

    /// Returns a reference to the `Kaula` struct if the dynamical tide is `KaulaTides`.
    pub fn kaula_get(&self) -> Option<&Kaula> {
        match &self.dynamical {
            DynamicalTide::KaulaTides(kaula) => Some(kaula),
            _ => None,
        }
    }

    /// Returns a mutable reference to the `Kaula` struct if the dynamical tide is `KaulaTides`.
    pub fn kaula_get_mut(&mut self) -> Option<&mut Kaula> {
        match &mut self.dynamical {
            DynamicalTide::KaulaTides(kaula) => Some(kaula),
            _ => None,
        }
    }
}

// `Tides` is read either from the keyword "Disabled" or from a map with the entries
// `equilibrium` and `dynamical`.
impl<'de> Deserialize<'de> for Tides {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields, default)]
        struct Entries {
            equilibrium: EquilibriumTide,
            dynamical: DynamicalTide,
        }
        impl Default for Entries {
            fn default() -> Self {
                Entries {
                    equilibrium: EquilibriumTide::Disabled,
                    dynamical: DynamicalTide::Disabled,
                }
            }
        }

        const LAYOUT: &str = "the tides of a body are given as `\"Disabled\"` or as `{\"equilibrium\": ..., \"dynamical\": ...}`; \
            the former `StellarConstantTimeLag`, `PlanetaryConstantTimeLag` and `KaulaTides` keys are now \
            `{\"equilibrium\": {\"ConstantTimeLag\": {...}}}`, `{\"dynamical\": \"FrequencyAveraged\"}` and \
            `{\"dynamical\": {\"KaulaTides\": {...}}}` (see docs/README.md)";

        struct TidesVisitor;

        impl<'de> Visitor<'de> for TidesVisitor {
            type Value = Tides;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("\"Disabled\" or a map with `equilibrium` and `dynamical`")
            }

            fn visit_str<E>(self, value: &str) -> std::result::Result<Tides, E>
            where
                E: de::Error,
            {
                if value == "Disabled" {
                    Ok(Tides::default())
                } else {
                    Err(E::custom(format!(
                        "unknown tides keyword `{value}`: {LAYOUT}"
                    )))
                }
            }

            fn visit_map<M>(self, map: M) -> std::result::Result<Tides, M::Error>
            where
                M: MapAccess<'de>,
            {
                let entries = Entries::deserialize(de::value::MapAccessDeserializer::new(map))
                    .map_err(|err| de::Error::custom(format!("{err}: {LAYOUT}")))?;
                Ok(Tides {
                    equilibrium: entries.equilibrium,
                    dynamical: entries.dynamical,
                })
            }
        }

        if deserializer.is_human_readable() {
            deserializer.deserialize_any(TidesVisitor)
        } else {
            // Binary formats (postcard) are not self-describing: the two entries are read
            // in the order they were written.
            let entries = Entries::deserialize(deserializer)?;
            Ok(Tides {
                equilibrium: entries.equilibrium,
                dynamical: entries.dynamical,
            })
        }
    }
}

#[cfg(test)]
pub mod tests;

// References:
// Benbakoura et al. 2019, https://doi.org/10.1051/0004-6361/201833314
// Bolmont & Mathis 2016, https://doi.org/10.1007/s10569-016-9690-3
// Kaula 1964, https://doi.org/10.1029/RG002i004p00661
