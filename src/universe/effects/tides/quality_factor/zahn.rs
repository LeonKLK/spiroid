use crate::constants::PI;
use crate::universe::particles::{Planet, Star};
use serde::{Deserialize, Serialize};

// Order unity parameters (Mustill & Villaver, 2012)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Zahn {
    f_prime: f64,
    c_f: f64,
    gamma_f: f64,
}

impl Default for Zahn {
    fn default() -> Self {
        // Values adopted from Mustill & Villaver (2012)
        Self {
            f_prime: 9. / 5., // multiplication factor for the dissipation
            c_f: 1.,
            gamma_f: 2., // exponent for the frequency dependence of the viscosity
        }
    }
}

impl Zahn {
    // Equilibrium tide with Zahn prescription
    // following the parametrisation (Mustill & Villaver, 2012, Eq. 1)
    pub(crate) fn tidal_quality(&self, star: &Star, planet: &Planet) -> f64 {
        let f2 = self.f_prime
            * min!(
                ((2. * PI) / (2. * planet.mean_motion * self.c_f * star.convective_turnover_time))
                    .powf(self.gamma_f),
                1.
            );
        let k2 = 1. / 27. / star.convective_turnover_time
            * (1. - star.radiative_mass / star.mass)
            * (star.mass + planet.mass)
            / star.mass
            / planet.mean_motion
            * (star.radius / planet.semi_major_axis).powi(3)
            * (2. * f2);

        3. / 2. / k2
    }
}

// References:
// Mustill & Villaver 2012, http://doi.org/10.1088/0004-637X/761/2/121
