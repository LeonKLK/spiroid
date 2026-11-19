use super::*;
use crate::Universe;
use crate::universe::effects::general_relativity::GeneralRelativityModel;
use crate::universe::effects::magnetism::{IsothermalWind, MagneticModel};
use crate::universe::effects::tides::constant_time_lag::tests::{
    test_planetary_constant_time_lag, test_stellar_constant_time_lag,
};
use crate::universe::effects::tides::kaula::tests::test_kaula;
use crate::universe::effects::tides::tests::test_stellar_tides;
use crate::universe::effects::tides::{ConstantTimeLag, DynamicalTide, EquilibriumTide, Tides};
use crate::universe::effects::wind::WindModel;
use crate::universe::particles::planet::tests::{
    test_planet, test_planet_kaula, test_planet_magnetic, test_planet_mercury,
};
use crate::universe::particles::star::tests_mode_ardestani::{test_star, test_star_evolving};
use crate::universe::tests::{DISK_IS_DISSIPATED, TEST_DISK_LIFETIME, TEST_TIME};
use crate::universe::{Particle, ParticleType};
use assert_approx_eq::assert_approx_eq;
use pretty_assertions::assert_eq;

#[test]
fn _derivatives_magnetic() {
    let star = test_star_evolving();
    let planet = test_planet_magnetic();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);

    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: Tides::default(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: Tides::default(),
            magnetism: MagneticModel::Wind(wind),
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };
    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -6.348994811695528e22;
    expected.central_body.convective_zone_angular_momentum = -9.407398756718968e21;
    expected.orbiting_body.semi_major_axis = -1.5957400309492905e43;

    assert_eq!(expected, result);
}

#[test]
fn _derivatives_tides() {
    let star = test_star_evolving();
    let planet = test_planet();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: Tides::default(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: test_stellar_tides(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -6.348994811695528e22;
    expected.central_body.convective_zone_angular_momentum = 5.72030250188552e23;
    expected.orbiting_body.semi_major_axis = -1.98486471648128e44;

    assert_eq!(expected, result);
}

#[test]
fn _derivatives_magnetic_tides() {
    let star = test_star_evolving();
    let planet = test_planet_magnetic();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);

    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: Tides::default(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: test_stellar_tides(),
            magnetism: MagneticModel::Wind(wind),
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.central_body.radiative_zone_angular_momentum = -6.348994811695528e22;
    expected.central_body.convective_zone_angular_momentum = 6.228617876044014e23;
    expected.orbiting_body.semi_major_axis = -2.144438719576209e44;

    assert_eq!(expected, result);
}

// Note: Tests are currently tailored for amd64 and may fail (with slight numerical discrepancies)
// on other platforms due to differences in cpu architecture.
#[test]
fn _derivatives_kaula() {
    let star = test_star();
    let planet = test_planet_kaula();

    let mut y = UniverseIntegral::default();
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);
    y.orbiting_body.spin = 8.062093352143078e-7;
    y.orbiting_body.eccentricity = 2.500000000179822e-5;
    y.orbiting_body.inclination = 0.34999207817863753;
    y.orbiting_body.longitude_ascending_node = 1.0465602799892118;
    y.orbiting_body.pericentre_omega = -0.11536773671287792;
    y.orbiting_body.spin_inclination = 0.31581363067032314;

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: Tides {
                dynamical: DynamicalTide::KaulaTides(test_kaula()),
                ..Default::default()
            },
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: Tides::default(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };

    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    let _ = force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let mut expected = UniverseIntegral::default();
    expected.orbiting_body.semi_major_axis = 3.0436830855775857e49;
    expected.orbiting_body.spin = -1.5432986377090881e-9;
    expected.orbiting_body.eccentricity = 5.129250165061513e-16;
    expected.orbiting_body.inclination = 0.0007011714730044864;
    expected.orbiting_body.longitude_ascending_node = -0.0018601514572935667;
    expected.orbiting_body.pericentre_omega = 5.876804234917257e-6;
    expected.orbiting_body.spin_inclination = 0.00035271712433232114;

    assert_eq!(expected, result);
}

// Constant time lag tide on both bodies: the star and the planet go through the same
// equations with the roles swapped, and each tide block of `force` writes its torque on
// its own deformed body (angular momentum of the convective zone for the star, spin for
// the planet).
fn constant_time_lag_universe(planet_spin_mode: PlanetarySpin) -> (Universe, UniverseIntegral) {
    let star = test_star_evolving();
    let mut planet = test_planet();
    planet.radius_of_gyration_2 = 0.254;
    planet.initialise();
    let spin = 3.0 * planet.mean_motion;
    let eccentricity: f64 = 0.3;

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    y.orbiting_body.semi_major_axis = planet.semi_major_axis.powf(6.5);
    y.orbiting_body.spin = spin;
    y.orbiting_body.eccentricity = eccentricity.powi(2);

    let universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: Tides {
                equilibrium: EquilibriumTide::ConstantTimeLag(ConstantTimeLag {
                    spin: planet_spin_mode,
                    ..test_planetary_constant_time_lag()
                }),
                dynamical: DynamicalTide::Disabled,
            },
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: Tides {
                equilibrium: EquilibriumTide::ConstantTimeLag(test_stellar_constant_time_lag()),
                dynamical: DynamicalTide::Disabled,
            },
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };
    (universe, y)
}

#[test]
fn _derivatives_constant_time_lag_star_and_planet() {
    let (mut universe, y) = constant_time_lag_universe(PlanetarySpin::Free);
    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let ParticleType::Star(star) = &universe.central_body.kind else {
        unreachable!()
    };
    let ParticleType::Planet(planet) = &universe.orbiting_body.kind else {
        unreachable!()
    };
    let stellar_ctl = universe.central_body.tides.constant_time_lag_get().unwrap();
    let planetary_ctl = universe
        .orbiting_body
        .tides
        .constant_time_lag_get()
        .unwrap();

    // Star: the torque of the stellar tide is added to the convective zone.
    let stellar_torque = stellar_ctl.tidal_torque(star, planet, planet);
    assert!(stellar_torque != 0.0);
    assert_eq!(
        star_convective_zone_angular_momentum_derivative(star, DISK_IS_DISSIPATED) + stellar_torque,
        result.central_body.convective_zone_angular_momentum
    );
    // The value recorded in the output is the same torque.
    assert_eq!(stellar_torque, star.tidal_torque_convective);

    // Planet: the torque of the planetary tide divided by the moment of inertia.
    let planetary_torque = planetary_ctl.tidal_torque(planet, star, planet);
    assert!(planetary_torque != 0.0);
    assert_eq!(
        planetary_torque / planet.moment_of_inertia,
        result.orbiting_body.spin
    );

    // The two torques are balanced by the orbit: d(L_orb)/dt = -(torque on the star) -
    // (torque on the planet), with L_orb = mu sqrt(G (M_star + M_p) a (1 - e^2)).
    let e2 = planet.eccentricity.powi(2);
    let orbital_angular_momentum = planet.reduced_mass
        * sqrt!(GRAVITATIONAL * (star.mass + planet.mass) * planet.semi_major_axis * (1. - e2));
    let semi_major_axis_rate =
        result.orbiting_body.semi_major_axis / (6.5 * planet.semi_major_axis.powf(5.5));
    let orbital_angular_momentum_rate = orbital_angular_momentum
        * (semi_major_axis_rate / (2. * planet.semi_major_axis)
            - result.orbiting_body.eccentricity / (2. * (1. - e2)));
    assert_approx_eq!(
        -(stellar_torque + planetary_torque) / orbital_angular_momentum_rate,
        1.0,
        1e-10
    );
}

#[test]
// With the planet spin held pseudo-synchronous the planetary torque is not integrated;
// the stellar tide is unaffected.
fn _derivatives_constant_time_lag_pseudo_synchronous() {
    let (mut universe, y) = constant_time_lag_universe(PlanetarySpin::PseudoSynchronous);
    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let ParticleType::Star(star) = &universe.central_body.kind else {
        unreachable!()
    };
    let ParticleType::Planet(planet) = &universe.orbiting_body.kind else {
        unreachable!()
    };
    let stellar_ctl = universe.central_body.tides.constant_time_lag_get().unwrap();

    assert_eq!(0.0, result.orbiting_body.spin);
    assert_eq!(
        ConstantTimeLag::pseudo_synchronous_spin(planet),
        planet.spin
    );
    assert_eq!(
        star_convective_zone_angular_momentum_derivative(star, DISK_IS_DISSIPATED)
            + stellar_ctl.tidal_torque(star, planet, planet),
        result.central_body.convective_zone_angular_momentum
    );
}

#[test]
// Two stellar tides together (constant time lag and frequency-averaged): their torques are
// summed first and the total is added once to the convective zone, so the derivative is
// exactly the non-tidal part plus the total tidal torque recorded in the output.
fn _derivatives_stellar_tides_collected() {
    let star = test_star_evolving();
    let planet = test_planet();

    let mut y = UniverseIntegral::default();
    y.central_body.radiative_zone_angular_momentum = star.radiative_zone_angular_momentum;
    y.central_body.convective_zone_angular_momentum = star.convective_zone_angular_momentum;
    // The frequency-averaged tide only acts when n < 2 spin (inertial waves): place the
    // planet far enough for that (0.3 AU; `test_planet` is inside this limit).
    let semi_major_axis: f64 = 4.5e10;
    y.orbiting_body.semi_major_axis = semi_major_axis.powf(6.5);

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: Tides::default(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Disabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: test_stellar_tides(),
            magnetism: MagneticModel::Disabled,
            wind: WindModel::Enabled,
            general_relativity: GeneralRelativityModel::Disabled,
        },
        time: TEST_TIME,
        disk_lifetime: TEST_DISK_LIFETIME,
        disk_is_dissipated: DISK_IS_DISSIPATED,
        derivatives: UniverseIntegral::default(),
    };
    universe.update(TEST_TIME, &y).unwrap();
    let mut result = UniverseIntegral::default();
    force(
        &universe.central_body,
        &universe.orbiting_body,
        universe.disk_is_dissipated,
        &mut result,
    )
    .unwrap();

    let ParticleType::Star(star) = &universe.central_body.kind else {
        unreachable!()
    };
    let ParticleType::Planet(planet) = &universe.orbiting_body.kind else {
        unreachable!()
    };
    let tides = &universe.central_body.tides;
    let quality_factor_torque = tides.stellar_quality_factor_torque(star, planet);
    let constant_time_lag_torque = tides
        .constant_time_lag_get()
        .unwrap()
        .tidal_torque(star, planet, planet);
    // Both tides are active in this configuration.
    assert!(planet.mean_motion < 2. * star.spin);
    assert!(quality_factor_torque != 0.0);
    assert!(constant_time_lag_torque != 0.0);
    assert_eq!(
        quality_factor_torque + constant_time_lag_torque,
        star.tidal_torque_convective
    );
    assert_eq!(
        star_convective_zone_angular_momentum_derivative(star, DISK_IS_DISSIPATED)
            + star.tidal_torque_convective,
        result.central_body.convective_zone_angular_momentum
    );
}

#[test]
fn _star_radiative_zone_angular_momentum_derivative() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_radiative_zone_angular_momentum_derivative(&star);
    let expected = -6.348994811695822e22;
    assert_eq!(expected, result);
}

#[test]
fn _star_convective_zone_angular_momentum_derivative() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_convective_zone_angular_momentum_derivative(&star, DISK_IS_DISSIPATED);
    let expected = -5.999284723052003e22;
    assert_eq!(expected, result);
}

#[test]
fn _planet_semi_major_axis_13_div_2_derivative() {
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let tides = test_stellar_tides();
    let mut wind = IsothermalWind::default();
    wind.footpoint_conductance = 7e4;
    let mut magnetism = MagneticModel::Wind(wind);
    let tidal_torque_convective = tides.stellar_tidal_torque(&star, &planet);
    let magnetic_torque = magnetism.magnetic_torque(&planet, &star);
    let wind_torque = WindModel::Enabled.wind_torque();

    star.update_wind_torque(wind_torque);
    star.update_tidal_torque(tidal_torque_convective);
    star.update_magnetic_torque(magnetic_torque);

    let result = planet_semi_major_axis_13_div_2_derivative_from_torque(
        &planet,
        &star,
        star.magnetic_torque + star.tidal_torque_convective,
    );
    let expected = -2.143941069104275e44;
    assert_eq!(expected, result);
}

#[test]
fn _kaula_planet_semi_major_axis_13_div_2_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);

    let tides = test_stellar_tides();
    let mut magnetism = MagneticModel::Wind(IsothermalWind::default());
    let tidal_torque_convective = tides.stellar_tidal_torque(&star, &planet);
    let magnetic_torque = magnetism.magnetic_torque(&planet, &star);
    let wind_torque = WindModel::Enabled.wind_torque();

    star.update_wind_torque(wind_torque);
    star.update_tidal_torque(tidal_torque_convective);
    star.update_magnetic_torque(magnetic_torque);

    let kaula = test_kaula();

    let result = kaula_planet_semi_major_axis_13_div_2_derivative(&planet, &star, &kaula);
    let expected = -2.3102674836008928e52;
    assert_eq!(expected, result);
}

#[test]
fn _planet_spin_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_derivative(&planet, &star, &kaula);
    let expected = 1.2501842317327892e-6;
    assert_eq!(expected, result);
}

#[test]
fn _planet_eccentricity_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_eccentricity_derivative(&planet, &star, &kaula);
    let expected = -2.371917363949444e-13;
    assert_eq!(expected, result);
}

#[test]
fn _planet_inclination_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.26702870058883815;
    assert_eq!(expected, result);
}

#[test]
fn _planet_longitude_ascending_node_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_longitude_ascending_node_derivative(&planet, &star, &kaula);
    let expected = 0.6637879161156786;
    assert_eq!(expected, result);
}

#[test]
fn _planet_argument_pericentre_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_argument_pericentre_derivative(&planet, &star, &kaula);
    let expected = -42049166.159453586;
    assert_eq!(expected, result);
}

#[test]
fn _kaula_star_argument_pericentre_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = kaula_star_argument_pericentre_derivative(&planet, &star, &kaula);
    let expected: f64 = -784585.0652653058;
    assert_approx_eq!(expected, result, abs!(expected) * 1e-12);
}

// In the coplanar case the stellar and planetary pericentre derivatives only differ by
// the torque prefactor of the deformed body: M_perturber^2 R_deformed^5.
#[test]
fn _kaula_star_argument_pericentre_derivative_symmetric() {
    let mut star = test_star();
    let mut planet = test_planet_kaula();
    // Coplanar: inclination and spin inclination to zero, keep spin, e and the angles.
    planet.refresh_orbital_elements(
        planet.spin,
        planet.eccentricity,
        0.0,
        planet.longitude_ascending_node,
        planet.pericentre_omega,
        0.0,
    );
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();

    let star_role = kaula_star_argument_pericentre_derivative(&planet, &star, &kaula);
    let planet_role = planet_argument_pericentre_derivative(&planet, &star, &kaula);
    assert_ne!(0.0, star_role);
    assert_approx_eq!(
        star_role * star.mass.powi(2) * planet.radius.powi(5),
        planet_role * planet.mass.powi(2) * star.radius.powi(5),
        abs!(planet_role * planet.mass.powi(2) * star.radius.powi(5)) * 1e-12
    );
}

#[test]
fn _planet_spin_axis_inclination_derivative() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_axis_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.1641163464603492;
    assert_eq!(expected, result);
}

#[test]
fn _general_relativity_pericentre_precession_rate() {
    let planet = test_planet_mercury();
    let result = general_relativity_pericentre_precession_rate(&planet);
    let expected = 4.725582990934671e-14;
    assert_eq!(expected, result);
}
