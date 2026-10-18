use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub enum WindModel {
    Disabled,
    // Matt et al. 2015 (Rossby number scaling), normalised by the star's `wind_torque_prefactor`.
    #[default]
    Enabled,
    // Bouvier et al. 1997 (as in Bolmont & Mathis 2016, Eq. 14, and Posidonius `wind.rs`):
    // dJ/dt = -k_wind * spin * min(spin, spin_saturation)^2 * (R / R_sun)^(1/2) * (M / M_sun)^(-1/2)
    // Hard switch at the saturation spin, no Rossby number dependence.
    Bouvier1997 {
        k_wind: f64,          // (kg.m2.s)
        spin_saturation: f64, // (rad.s-1)
    },
}

impl WindModel {
    // Returns `true` if a stellar wind torque is to be computed (any model but `Disabled`).
    pub(crate) fn wind_torque(&self) -> bool {
        *self != Self::Disabled
    }
}
