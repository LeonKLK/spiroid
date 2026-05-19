use crate::constants::{SECONDS_IN_YEAR, SOLAR_LUMINOSITY, SOLAR_MASS, SOLAR_RADIUS};
use derive_more::{Add, Mul};
use serde::{Deserialize, Serialize};

// Interpolation values deserialized from user provided CSV.
// Source of data is typically from STAREVOL or MESA stellar models.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone, Copy, Add, Mul)]
pub struct StarCsv {
    pub(crate) age: f64,                          // (s)
    pub(crate) radius: f64,                       // (m)
    pub(crate) mass: f64,                         // (kg)
    pub(crate) convective_radius: f64,            // (m)
    pub(crate) radiative_mass: f64,               // (kg)
    pub(crate) radiative_moment_of_inertia: f64,  // (kg.m2)
    pub(crate) convective_moment_of_inertia: f64, // (kg.m2)
    pub(crate) luminosity: f64,                   // (J.s-1)

    // Only provided by MESA data files.
    #[serde(default)]
    pub(crate) convective_turnover_time: f64, // (s)
    #[serde(default)]
    pub(crate) mass_loss_rate: f64, // (kg.s-1)

    // Calculated internally, not included in the CSV.
    #[serde(skip)]
    pub(crate) convective_moment_of_inertia_derivative: f64,
    #[serde(skip)]
    pub(crate) radiative_mass_derivative: f64,
    #[serde(skip)]
    pub(crate) radius_derivative: f64,
}

impl StarCsv {
    pub fn initialise(stars: &mut [Self]) {
        stars.iter_mut().for_each(Self::convert_units);
        Self::compute_derivatives(stars);
    }

    // Return the ages, which are used as the index to interpolate remaining stellar values, based on time.
    pub fn ages(stars: &[Self]) -> Vec<f64> {
        stars
            .iter()
            .map(|starcsv| starcsv.age)
            .collect::<Vec<f64>>()
    }

    // Initialise the input values with unit conversion.
    fn convert_units(&mut self) {
        self.age *= SECONDS_IN_YEAR;
        self.radius *= SOLAR_RADIUS;
        self.mass *= SOLAR_MASS;
        self.luminosity *= SOLAR_LUMINOSITY;
        self.convective_radius *= SOLAR_RADIUS;
        self.radiative_mass *= SOLAR_MASS;
        self.radiative_moment_of_inertia *= self.mass * self.radius.powi(2);
        // convective_moment_of_inertia is now the moment of inertia of the whole solid body star
        self.convective_moment_of_inertia = 0.07 * self.mass * self.radius.powi(2);
        self.mass_loss_rate *= SOLAR_MASS / SECONDS_IN_YEAR;
    }

    // Calcultes the radiative_mass_derivative and convective_moment_of_inertia_derivative for each record.
    fn compute_derivatives(stars: &mut [Self]) {
        let stars_len = stars.len();

        // Derivative is zero for first and last timesteps.
        stars[0].radiative_mass_derivative = 0.;
        stars[0].convective_moment_of_inertia_derivative = 0.;
        stars[0].radius_derivative = 0.;
        stars[stars_len - 1].radiative_mass_derivative = 0.;
        stars[stars_len - 1].convective_moment_of_inertia_derivative = 0.;
        stars[stars_len - 1].radius_derivative = 0.;

        for i in 1..stars_len - 1 {
            // Unpack values of the star at three consecutive timesteps to compute the derivatives.
            let [prev, curr, next] = &mut stars[i - 1..=i + 1] else {
                unreachable!()
            };
            curr.radiative_mass_derivative =
                (next.radiative_mass - prev.radiative_mass) / (next.age - prev.age);
            curr.convective_moment_of_inertia_derivative = (next.convective_moment_of_inertia
                - prev.convective_moment_of_inertia)
                / (next.age - prev.age);
            curr.radius_derivative = (next.radius - prev.radius) / (next.age - prev.age);
        }
    }
}
