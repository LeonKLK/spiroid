use super::*;

use crate::universe::particles::planet::tests::test_planet;
use crate::universe::particles::star::tests_mode_ardestani::test_star;

use pretty_assertions::assert_eq;

fn assert_close(expected: f64, result: f64, relative_tolerance: f64) {
    assert!(
        abs!(expected - result) <= relative_tolerance * abs!(expected),
        "expected {expected}, got {result}"
    );
}

#[test]
// Benbakoura et al. (2019) Eq. 19: the torque is inversely proportional to Q' and takes the
// sign opposite to the tidal frequency 2 (omega_star - n).
fn _stellar_tidal_torque() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    // The test star spins slower than the orbit: it is spun up.
    assert!(star.tidal_frequency < 0.0);

    let tidal_quality = 1e6;
    let expected = (9. / 4.) * GRAVITATIONAL * planet.mass.powi(2) * star.radius.powi(5)
        / planet.semi_major_axis.powi(6)
        / tidal_quality;
    let result = stellar_tidal_torque(tidal_quality, &star, &planet);
    assert_close(expected, result, 1e-12);
    assert_close(
        result / 10.,
        stellar_tidal_torque(10. * tidal_quality, &star, &planet),
        1e-12,
    );

    // A star spinning faster than the orbit is spun down.
    star.spin = 2.0 * planet.mean_motion;
    star.refresh_tidal_frequency(&planet);
    assert_close(
        -expected,
        stellar_tidal_torque(tidal_quality, &star, &planet),
        1e-12,
    );

    // No torque at corotation.
    star.tidal_frequency = 0.0;
    assert_eq!(0.0, stellar_tidal_torque(tidal_quality, &star, &planet));
}

#[test]
// Mustill & Villaver (2012) Eq. 1 with their order unity parameters.
fn _zahn_tidal_quality() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);

    let zahn = Zahn::default();
    let result = zahn.tidal_quality(&star, &planet);
    assert!(result.is_finite() && result > 0.0);

    // The dissipation is proportional to the multiplication factor f'.
    let json = r#"{"f_prime": 3.6, "c_f": 1.0, "gamma_f": 2.0}"#;
    let doubled: Zahn = serde_json::from_str(json).unwrap();
    assert_close(result / 2., doubled.tidal_quality(&star, &planet), 1e-12);
}

#[test]
// Q'_dyn = 3 / (2 <D>_omega), <D>_omega = dissipation * omega_star^2 (Bolmont & Mathis 2016
// Eq. 4), while inertial waves can be excited (n < 2 omega_star); no dissipation otherwise.
fn _frequency_averaged_tidal_quality() {
    let mut star = test_star();
    let planet = test_planet();
    star.dynamical_tide_dissipation = 3037.223911926055;

    // Inertial waves are excited: tidal frequency above -2 omega_star.
    star.spin = planet.mean_motion;
    star.refresh_tidal_frequency(&planet);
    let expected = 3. / (2. * star.dynamical_tide_dissipation * star.spin.powi(2));
    assert_close(expected, frequency_averaged::tidal_quality(&star), 1e-12);

    // Outside the inertial range the tide is switched off (Q' goes to infinity).
    star.spin = 0.1 * planet.mean_motion;
    star.refresh_tidal_frequency(&planet);
    assert!(frequency_averaged::tidal_quality(&star) > 1e6 * expected);
}
