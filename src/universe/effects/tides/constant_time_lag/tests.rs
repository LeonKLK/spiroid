use super::*;

// Stellar constant time lag tide
use crate::universe::particles::planet::tests::test_planet;
use crate::universe::particles::star::tests_mode_ardestani::test_star;

use pretty_assertions::assert_eq;

pub fn test_constant_time_lag() -> StellarConstantTimeLag {
    StellarConstantTimeLag {
        equilibrium: Equilibrium::SigmaBarStar(1e-6),
        inertial: Inertial::FrequencyAveraged,
    }
}

#[test]
// This function is only called if tides are enabled.
fn _tidal_torque() {
    let expected = 6.325284391272144e23;
    let mut star = test_star();
    let planet = test_planet();
    let tidal_model = test_constant_time_lag();
    star.refresh_tidal_frequency(&planet);
    let result = tidal_model.stellar_tidal_torque(&star, &planet);
    assert_eq!(expected, result);
}

// Planetary constant time lag tide

use crate::constants::GRAVITATIONAL;

pub fn test_planetary_constant_time_lag() -> PlanetaryConstantTimeLag {
    PlanetaryConstantTimeLag {
        // Jupiter, Leconte et al. (2010)
        k2_delta_t: 0.025,
        spin: PlanetarySpin::Free,
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

#[test]
fn _planetary_eccentricity_functions_circular() {
    assert_eq!(1.0, PlanetaryConstantTimeLag::big_n(0.0));
    assert_eq!(1.0, PlanetaryConstantTimeLag::big_n_a(0.0));
    assert_eq!(1.0, PlanetaryConstantTimeLag::big_omega_e(0.0));
    assert_eq!(1.0, PlanetaryConstantTimeLag::big_n_e(0.0));
    assert_eq!(1.0, PlanetaryConstantTimeLag::big_omega(0.0));
}

#[test]
// Hut (1981) polynomials evaluated at e = 0.3 (Leconte et al. 2010 Eqs. 3, 4, 7, 8, 11).
fn _planetary_eccentricity_functions() {
    let e = 0.3;
    assert_close(
        3.0302579446153191e+00,
        PlanetaryConstantTimeLag::big_n(e),
        1e-14,
    );
    assert_close(
        5.3993220183297863e+00,
        PlanetaryConstantTimeLag::big_n_a(e),
        1e-14,
    );
    assert_close(
        1.8204390734037521e+00,
        PlanetaryConstantTimeLag::big_omega_e(e),
        1e-14,
    );
    assert_close(
        2.4971619665511264e+00,
        PlanetaryConstantTimeLag::big_n_e(e),
        1e-14,
    );
    assert_close(
        1.9460541987451936e+00,
        PlanetaryConstantTimeLag::big_omega(e),
        1e-14,
    );
}

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
    assert_close(expected, ctl.planetary_tidal_torque(&planet, &star), 1e-12);
}

#[test]
// Circular orbit: the angular momentum lost by the planet spin goes into the orbit,
// C_p d(omega_p)/dt + dL_orb/dt = 0 with L_orb = mu sqrt(G (M_star + M_p) a).
fn _planetary_angular_momentum_conservation_circular() {
    let star = test_star();
    let planet = test_planet_ctl(3.0 * test_planet().mean_motion, 0.0);
    let ctl = test_planetary_constant_time_lag();
    let spin_angular_momentum_rate = planet.moment_of_inertia * ctl.spin_derivative(&planet, &star);
    // d(a^13/2)/dt = 13/2 a^11/2 da/dt
    let semi_major_axis_rate = ctl.semi_major_axis_13_div_2_derivative(&planet, &star)
        / (6.5 * planet.semi_major_axis.powf(5.5));
    let orbital_angular_momentum = planet.reduced_mass
        * sqrt!(GRAVITATIONAL * (star.mass + planet.mass) * planet.semi_major_axis);
    let orbital_angular_momentum_rate =
        orbital_angular_momentum / (2.0 * planet.semi_major_axis) * semi_major_axis_rate;
    assert_close(
        -spin_angular_momentum_rate,
        orbital_angular_momentum_rate,
        1e-10,
    );
}

#[test]
// At the pseudo-synchronous spin the spin derivative vanishes (Leconte et al. 2010 Eq. 12).
fn _planetary_pseudo_synchronous_spin() {
    let star = test_star();
    let mut planet = test_planet_ctl(0.0, 0.3);
    planet.spin = PlanetaryConstantTimeLag::pseudo_synchronous_spin(&planet);
    assert!(planet.spin > planet.mean_motion);
    let ctl = test_planetary_constant_time_lag();
    let scale = abs!(ctl.prefactor(&planet, &star) / planet.mean_motion / planet.moment_of_inertia);
    assert!(abs!(ctl.spin_derivative(&planet, &star)) < 1e-12 * scale);
    // For e = 0 the pseudo-synchronous spin is the mean motion.
    let circular = test_planet_ctl(0.0, 0.0);
    assert_eq!(
        circular.mean_motion,
        PlanetaryConstantTimeLag::pseudo_synchronous_spin(&circular)
    );
}

#[test]
// A pseudo-synchronous planet on an eccentric orbit circularises and migrates inward.
fn _planetary_tide_circularises() {
    let star = test_star();
    let mut planet = test_planet_ctl(0.0, 0.3);
    planet.spin = PlanetaryConstantTimeLag::pseudo_synchronous_spin(&planet);
    let ctl = PlanetaryConstantTimeLag {
        spin: PlanetarySpin::PseudoSynchronous,
        ..test_planetary_constant_time_lag()
    };
    assert!(ctl.eccentricity_squared_derivative(&planet, &star) < 0.0);
    assert!(ctl.semi_major_axis_13_div_2_derivative(&planet, &star) < 0.0);
    // The spin is not integrated in this mode.
    assert_eq!(0.0, ctl.spin_derivative(&planet, &star));
    // No eccentricity change on a circular orbit.
    let circular = test_planet_ctl(0.0, 0.0);
    assert_eq!(0.0, ctl.eccentricity_squared_derivative(&circular, &star));
}
