use crate::universe::particles::Star;

// Dynamical tide
// Equivalent Q' factor, tidal quality factor
// This is the inverse of Eq. 4 of Bolmont & Mathis 2016,
// which give the tidal quality factor Q' as a function of <D>w
// (Eq. 1 of Mathis 2015, or Bolmont & Mathis 2016)
// But here, keep in mind that the spin was removed out of <D>w to be then multiplied here
// Q' = 3 / ( 2 * <D>w ) | Bolmont & Mathis
// Q' = 3 / ( 2 * <D>w' * spin^2 ) | here
pub(crate) fn tidal_quality(star: &Star) -> f64 {
    // Epsilon to ensure a smooth transition equilibrium / dynamical tide
    let epsilon_step = 1e-6;
    // dynamical_tide_quality_factors
    3. / (2. * star.dynamical_tide_dissipation * star.spin.powi(2)) * 2.
        / (1. + tanh!((star.tidal_frequency + 2. * star.spin) / epsilon_step))
}

// References:
// Mathis 2015, https://doi.org/10.1051/0004-6361/201526472
// Bolmont & Mathis 2016, https://doi.org/10.1007/s10569-016-9690-3
