use super::*;

use crate::constants::{GRAVITATIONAL, SOLAR_MASS, SOLAR_RADIUS};
use crate::universe::effects::tides::quality_factor;
use crate::universe::particles::Star;
use crate::universe::particles::planet::tests::test_planet;
use crate::universe::particles::star::tests_mode_ardestani::test_star;
use crate::universe::physics::planet_semi_major_axis_13_div_2_derivative_from_torque;

use pretty_assertions::assert_eq;

// Stellar constant time lag tide, given as the normalised dissipation of Hansen (2010).
pub fn test_stellar_constant_time_lag() -> ConstantTimeLag {
    ConstantTimeLag {
        sigma_bar_star: 1e-6,
        ..Default::default()
    }
}

// Planetary constant time lag tide.
pub fn test_planetary_constant_time_lag() -> ConstantTimeLag {
    ConstantTimeLag {
        // Jupiter, Leconte et al. (2010)
        k2_delta_t: 0.025,
        spin: PlanetarySpin::Free,
        ..Default::default()
    }
}

// `test_planet` with a spin, an eccentricity and a moment of inertia (coplanar).
fn test_planet_ctl(spin: f64, eccentricity: f64) -> Planet {
    let mut planet = test_planet();
    planet.radius_of_gyration_2 = 0.254;
    planet.initialise();
    planet.refresh_orbital_elements(spin, eccentricity, 0.0, 0.0, 0.0, 0.0);
    planet
}

fn assert_close(expected: f64, result: f64, relative_tolerance: f64) {
    assert!(
        abs!(expected - result) <= relative_tolerance * abs!(expected),
        "expected {expected}, got {result}"
    );
}

// Rate of change of the orbital angular momentum L = mu sqrt(G (M_star + M_p) a (1 - e^2))
// from the derivatives of the integrated variables a^13/2 and e^2:
// dL/dt = L [ (da/dt) / (2 a) - (de^2/dt) / (2 (1 - e^2)) ].
fn orbital_angular_momentum_rate(
    star: &Star,
    planet: &Planet,
    semi_major_axis_13_div_2_derivative: f64,
    eccentricity_squared_derivative: f64,
) -> f64 {
    let e2 = planet.eccentricity.powi(2);
    let orbital_angular_momentum = planet.reduced_mass
        * sqrt!(GRAVITATIONAL * (star.mass + planet.mass) * planet.semi_major_axis * (1. - e2));
    // d(a^13/2)/dt = 13/2 a^11/2 da/dt
    let semi_major_axis_rate =
        semi_major_axis_13_div_2_derivative / (6.5 * planet.semi_major_axis.powf(5.5));
    orbital_angular_momentum
        * (semi_major_axis_rate / (2.0 * planet.semi_major_axis)
            - eccentricity_squared_derivative / (2.0 * (1. - e2)))
}

#[test]
fn _eccentricity_functions_circular() {
    assert_eq!(1.0, ConstantTimeLag::big_n(0.0));
    assert_eq!(1.0, ConstantTimeLag::big_n_a(0.0));
    assert_eq!(1.0, ConstantTimeLag::big_omega_e(0.0));
    assert_eq!(1.0, ConstantTimeLag::big_n_e(0.0));
    assert_eq!(1.0, ConstantTimeLag::big_omega(0.0));
}

#[test]
// Hut (1981) polynomials evaluated at e = 0.3 (Leconte et al. 2010 Eqs. 3, 4, 7, 8, 11).
fn _eccentricity_functions() {
    let e = 0.3;
    assert_close(3.0302579446153191e+00, ConstantTimeLag::big_n(e), 1e-14);
    assert_close(5.3993220183297863e+00, ConstantTimeLag::big_n_a(e), 1e-14);
    assert_close(
        1.8204390734037521e+00,
        ConstantTimeLag::big_omega_e(e),
        1e-14,
    );
    assert_close(2.4971619665511264e+00, ConstantTimeLag::big_n_e(e), 1e-14);
    assert_close(1.9460541987451936e+00, ConstantTimeLag::big_omega(e), 1e-14);
}

// Time lag of the deformed body

#[test]
// A direct `k2_delta_t` is used as given, for either body.
fn _k2_delta_t_direct() {
    let ctl = test_planetary_constant_time_lag();
    assert_eq!(0.025, ctl.k2_delta_t(&test_planet()));
    assert_eq!(0.025, ctl.k2_delta_t(&test_star()));
}

#[test]
// Bolmont & Mathis (2016) Eq. 8: k2 dt = 3/2 sigma_bar_star sigma_0 R^5 / G, i.e. about
// 2389 s * sigma_bar_star for a star of one solar radius, growing as R^5.
fn _k2_delta_t_from_sigma_bar_star() {
    let ctl = test_stellar_constant_time_lag();
    let mut star = test_star();

    star.radius = SOLAR_RADIUS;
    let sigma_0 = sqrt!(GRAVITATIONAL / (SOLAR_MASS * SOLAR_RADIUS.powi(7)));
    let expected = 1.5 * 1e-6 * sigma_0 * SOLAR_RADIUS.powi(5) / GRAVITATIONAL;
    assert_close(expected, ctl.k2_delta_t(&star), 1e-14);
    assert_close(2.389e-3, ctl.k2_delta_t(&star), 1e-3);

    star.radius = 2.0 * SOLAR_RADIUS;
    assert_close(32.0 * expected, ctl.k2_delta_t(&star), 1e-14);
}

// Planetary tide: planet deformed, star perturber.

#[test]
// For e = 0 the torque on the planet is Hut's -3 k2 dt G M_star^2 R_p^5 (omega_p - n) / a^6.
fn _planetary_tidal_torque_circular() {
    let star = test_star();
    let planet = test_planet_ctl(3.0 * test_planet().mean_motion, 0.0);
    let ctl = test_planetary_constant_time_lag();
    let expected = -3.0
        * ctl.k2_delta_t
        * GRAVITATIONAL
        * star.mass.powi(2)
        * planet.radius.powi(5)
        * (planet.spin - planet.mean_motion)
        / planet.semi_major_axis.powi(6);
    assert_close(expected, ctl.tidal_torque(&planet, &star, &planet), 1e-12);
}

#[test]
// The angular momentum lost by the planet spin goes into the orbit, at any eccentricity:
// C_p d(omega_p)/dt + dL_orb/dt = 0.
fn _planetary_angular_momentum_conservation() {
    let star = test_star();
    let ctl = test_planetary_constant_time_lag();
    for eccentricity in [0.0, 0.3, 0.7] {
        let planet = test_planet_ctl(3.0 * test_planet().mean_motion, eccentricity);
        let spin_angular_momentum_rate = ctl.tidal_torque(&planet, &star, &planet);
        let orbital_angular_momentum_rate = orbital_angular_momentum_rate(
            &star,
            &planet,
            ctl.semi_major_axis_13_div_2_derivative(&planet, &star, &planet),
            ctl.eccentricity_squared_derivative(&planet, &star, &planet),
        );
        assert_close(
            -spin_angular_momentum_rate,
            orbital_angular_momentum_rate,
            1e-10,
        );
    }
}

#[test]
// At the pseudo-synchronous spin the spin derivative vanishes (Leconte et al. 2010 Eq. 12).
fn _planetary_pseudo_synchronous_spin() {
    let star = test_star();
    let mut planet = test_planet_ctl(0.0, 0.3);
    planet.spin = ConstantTimeLag::pseudo_synchronous_spin(&planet);
    assert!(planet.spin > planet.mean_motion);
    let ctl = test_planetary_constant_time_lag();
    let scale = abs!(
        ctl.prefactor(&planet, &star, &planet) / planet.mean_motion / planet.moment_of_inertia
    );
    let spin_derivative = ctl.tidal_torque(&planet, &star, &planet) / planet.moment_of_inertia;
    assert!(abs!(spin_derivative) < 1e-12 * scale);
    // For e = 0 the pseudo-synchronous spin is the mean motion.
    let circular = test_planet_ctl(0.0, 0.0);
    assert_eq!(
        circular.mean_motion,
        ConstantTimeLag::pseudo_synchronous_spin(&circular)
    );
}

#[test]
// A pseudo-synchronous planet on an eccentric orbit circularises and migrates inward.
fn _planetary_tide_circularises() {
    let star = test_star();
    let mut planet = test_planet_ctl(0.0, 0.3);
    planet.spin = ConstantTimeLag::pseudo_synchronous_spin(&planet);
    let ctl = ConstantTimeLag {
        spin: PlanetarySpin::PseudoSynchronous,
        ..test_planetary_constant_time_lag()
    };
    assert!(ctl.eccentricity_squared_derivative(&planet, &star, &planet) < 0.0);
    assert!(ctl.semi_major_axis_13_div_2_derivative(&planet, &star, &planet) < 0.0);
    // No eccentricity change on a circular orbit.
    let circular = test_planet_ctl(0.0, 0.0);
    assert_eq!(
        0.0,
        ctl.eccentricity_squared_derivative(&circular, &star, &circular)
    );
}

// Stellar tide: star deformed, planet perturber. Same equations, roles swapped.

#[test]
// For e = 0 the torque on the star is Hut's -3 k2 dt G M_p^2 R_star^5 (omega_star - n) / a^6.
fn _stellar_tidal_torque_circular() {
    let star = test_star();
    let planet = test_planet();
    let ctl = test_stellar_constant_time_lag();
    let expected = -3.0
        * ctl.k2_delta_t(&star)
        * GRAVITATIONAL
        * planet.mass.powi(2)
        * star.radius.powi(5)
        * (star.spin - planet.mean_motion)
        / planet.semi_major_axis.powi(6);
    assert_close(expected, ctl.tidal_torque(&star, &planet, &planet), 1e-12);
}

#[test]
// For e = 0 the constant time lag torque on the star is the torque of the tidal quality
// factor framework (Benbakoura et al. 2019 Eq. 19) with the equilibrium tide quality factor
// Q'_eq = G / (sigma_star R^5 |omega_tide|), omega_tide = 2 (omega_star - n): the model that
// `sigma_bar_star` selected before the equilibrium tide moved to the constant time lag
// equations. The two only differ by the smoothing of the sign of the tidal frequency.
fn _stellar_tidal_torque_circular_is_quality_factor_torque() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let ctl = test_stellar_constant_time_lag();

    let sigma_star =
        ctl.sigma_bar_star * sqrt!(GRAVITATIONAL / (SOLAR_MASS * SOLAR_RADIUS.powi(7)));
    let tidal_quality =
        GRAVITATIONAL / (sigma_star * star.radius.powi(5) * abs!(star.tidal_frequency));
    let expected = quality_factor::stellar_tidal_torque(tidal_quality, &star, &planet);

    assert_close(expected, ctl.tidal_torque(&star, &planet, &planet), 1e-12);
}

#[test]
// For e = 0 the semi-major axis derivative of the constant time lag stellar tide is the one
// derived from its torque by angular momentum conservation (Benbakoura et al. 2019 Eq. 3).
fn _stellar_semi_major_axis_circular_is_torque_route() {
    let star = test_star();
    let planet = test_planet();
    let ctl = test_stellar_constant_time_lag();
    assert_eq!(0.0, star.evolved_change_semi_major_axis);
    let expected = planet_semi_major_axis_13_div_2_derivative_from_torque(
        &planet,
        &star,
        ctl.tidal_torque(&star, &planet, &planet),
    );
    assert_close(
        expected,
        ctl.semi_major_axis_13_div_2_derivative(&star, &planet, &planet),
        1e-12,
    );
    assert_eq!(
        0.0,
        ctl.eccentricity_squared_derivative(&star, &planet, &planet)
    );
}

#[test]
// The angular momentum taken from the stellar spin goes into the orbit, at any eccentricity:
// torque on the star + dL_orb/dt = 0.
fn _stellar_angular_momentum_conservation() {
    let star = test_star();
    let ctl = test_stellar_constant_time_lag();
    for eccentricity in [0.0, 0.3, 0.7] {
        let planet = test_planet_ctl(0.0, eccentricity);
        let orbital_angular_momentum_rate = orbital_angular_momentum_rate(
            &star,
            &planet,
            ctl.semi_major_axis_13_div_2_derivative(&star, &planet, &planet),
            ctl.eccentricity_squared_derivative(&star, &planet, &planet),
        );
        assert_close(
            -ctl.tidal_torque(&star, &planet, &planet),
            orbital_angular_momentum_rate,
            1e-10,
        );
    }
}

#[test]
// A star spinning slower than the orbit takes angular momentum from it: the planet migrates
// inward and the orbit circularises. A star spinning faster than 18/11 n (to lowest order in
// e) pushes the planet outward and raises the eccentricity (Leconte et al. 2010 Eq. 6).
fn _stellar_tide_direction() {
    let ctl = test_stellar_constant_time_lag();
    let planet = test_planet_ctl(0.0, 0.05);

    let mut slow = test_star();
    slow.spin = 0.5 * planet.mean_motion;
    assert!(ctl.tidal_torque(&slow, &planet, &planet) > 0.0);
    assert!(ctl.semi_major_axis_13_div_2_derivative(&slow, &planet, &planet) < 0.0);
    assert!(ctl.eccentricity_squared_derivative(&slow, &planet, &planet) < 0.0);

    let mut fast = test_star();
    fast.spin = 2.0 * planet.mean_motion;
    assert!(ctl.tidal_torque(&fast, &planet, &planet) < 0.0);
    assert!(ctl.semi_major_axis_13_div_2_derivative(&fast, &planet, &planet) > 0.0);
    assert!(ctl.eccentricity_squared_derivative(&fast, &planet, &planet) > 0.0);
}
