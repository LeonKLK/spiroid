use super::*;

use crate::universe::effects::tides::constant_time_lag::tests::{
    test_planetary_constant_time_lag, test_stellar_constant_time_lag,
};
use crate::universe::effects::{GeneralRelativityModel, MagneticModel, WindModel};
use crate::universe::particles::planet::tests::test_planet;
use crate::universe::particles::star::tests_mode_ardestani::test_star;
use crate::universe::particles::{Particle, ParticleType};

use pretty_assertions::assert_eq;

fn assert_close(expected: f64, result: f64, relative_tolerance: f64) {
    assert!(
        abs!(expected - result) <= relative_tolerance * abs!(expected),
        "expected {expected}, got {result}"
    );
}

// Constant time lag equilibrium tide and frequency-averaged dynamical tide on the star.
pub fn test_stellar_tides() -> Tides {
    Tides {
        equilibrium: EquilibriumTide::ConstantTimeLag(test_stellar_constant_time_lag()),
        dynamical: DynamicalTide::FrequencyAveraged,
    }
}

fn test_particle(kind: ParticleType, tides: Tides) -> Particle {
    Particle {
        kind,
        tides,
        magnetism: MagneticModel::Disabled,
        wind: WindModel::Disabled,
        general_relativity: GeneralRelativityModel::Disabled,
    }
}

fn initialise_error(mut particle: Particle) -> String {
    particle.initialise(0.0).unwrap_err().to_string()
}

// Torque on the star

#[test]
// No tidal torque if tides are disabled.
fn _stellar_tidal_torque_disabled() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let tides = Tides::default();
    assert_eq!(0.0, tides.stellar_quality_factor_torque(&star, &planet));
    assert_eq!(0.0, tides.stellar_tidal_torque(&star, &planet));
}

#[test]
// The torques of the equilibrium and dynamical tides add up (Benbakoura et al. 2019 Eq. 20).
// Only the quality factor tides go through the circular torque route of physics.rs.
fn _stellar_tidal_torque_is_sum_of_tides() {
    let mut star = test_star();
    let mut planet = test_planet();
    // Inside the range where inertial waves are excited (n < 2 omega_star).
    star.spin = 0.8 * planet.mean_motion;
    star.dynamical_tide_dissipation = 3037.223911926055;
    star.refresh_tidal_frequency(&planet);
    planet.refresh_orbital_elements(0.0, 0.2, 0.0, 0.0, 0.0, 0.0);

    let equilibrium = Tides {
        equilibrium: EquilibriumTide::ConstantTimeLag(test_stellar_constant_time_lag()),
        dynamical: DynamicalTide::Disabled,
    };
    let dynamical = Tides {
        equilibrium: EquilibriumTide::Disabled,
        dynamical: DynamicalTide::FrequencyAveraged,
    };
    let both = test_stellar_tides();

    assert!(dynamical.stellar_tidal_torque(&star, &planet) > 0.0);
    assert_eq!(
        0.0,
        equilibrium.stellar_quality_factor_torque(&star, &planet)
    );
    assert_close(
        test_stellar_constant_time_lag().tidal_torque(&star, &planet, &planet),
        equilibrium.stellar_tidal_torque(&star, &planet),
        1e-14,
    );
    assert_close(
        dynamical.stellar_tidal_torque(&star, &planet),
        both.stellar_quality_factor_torque(&star, &planet),
        1e-14,
    );
    assert_close(
        equilibrium.stellar_tidal_torque(&star, &planet)
            + dynamical.stellar_tidal_torque(&star, &planet),
        both.stellar_tidal_torque(&star, &planet),
        1e-14,
    );
}

#[test]
// Circular orbit: the torque of the former `StellarConstantTimeLag` model (equilibrium tide
// `SigmaBarStar(1e-6)` as a tidal quality factor, plus `FrequencyAveraged`) is recovered.
// 6.325284391272144e23 is the value that model gave for this star and planet.
fn _stellar_tidal_torque_matches_former_quality_factor_model() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = test_stellar_tides().stellar_tidal_torque(&star, &planet);
    assert_close(6.325284391272144e23, result, 1e-4);
}

#[test]
// The Zahn equilibrium tide and the frequency-averaged dynamical tide share one quality
// factor: 1 / Q' = 1 / Q'_eq + 1 / Q'_dyn.
fn _stellar_quality_factor_torque_zahn_and_frequency_averaged() {
    let mut star = test_star();
    let planet = test_planet();
    star.spin = 0.8 * planet.mean_motion;
    star.dynamical_tide_dissipation = 3037.223911926055;
    star.refresh_tidal_frequency(&planet);

    let zahn = Tides {
        equilibrium: EquilibriumTide::Zahn(Zahn::default()),
        dynamical: DynamicalTide::Disabled,
    };
    let dynamical = Tides {
        equilibrium: EquilibriumTide::Disabled,
        dynamical: DynamicalTide::FrequencyAveraged,
    };
    let both = Tides {
        equilibrium: EquilibriumTide::Zahn(Zahn::default()),
        dynamical: DynamicalTide::FrequencyAveraged,
    };
    assert_close(
        zahn.stellar_quality_factor_torque(&star, &planet),
        zahn.stellar_tidal_torque(&star, &planet),
        1e-14,
    );
    assert_close(
        zahn.stellar_tidal_torque(&star, &planet) + dynamical.stellar_tidal_torque(&star, &planet),
        both.stellar_tidal_torque(&star, &planet),
        1e-12,
    );
}

// Configuration file layout

#[test]
fn _deserialise_disabled() {
    let tides: Tides = serde_json::from_str(r#""Disabled""#).unwrap();
    assert_eq!(Tides::default(), tides);
    let tides: Tides = serde_json::from_str("{}").unwrap();
    assert_eq!(Tides::default(), tides);
}

#[test]
fn _deserialise_stellar_tides() {
    let json = r#"{
        "equilibrium": {"ConstantTimeLag": {"sigma_bar_star": 1e-06}},
        "dynamical": "FrequencyAveraged"
    }"#;
    let tides: Tides = serde_json::from_str(json).unwrap();
    assert_eq!(test_stellar_tides(), tides);

    // A missing entry is disabled.
    let tides: Tides = serde_json::from_str(r#"{"dynamical": "FrequencyAveraged"}"#).unwrap();
    assert_eq!(EquilibriumTide::Disabled, tides.equilibrium);
    assert_eq!(DynamicalTide::FrequencyAveraged, tides.dynamical);

    let tides: Tides = serde_json::from_str(
        r#"{"equilibrium": {"Zahn": {"f_prime": 1.8, "c_f": 1.0, "gamma_f": 2.0}}}"#,
    )
    .unwrap();
    assert_eq!(EquilibriumTide::Zahn(Zahn::default()), tides.equilibrium);
    assert_eq!(DynamicalTide::Disabled, tides.dynamical);
}

#[test]
fn _deserialise_planetary_constant_time_lag() {
    let json = r#"{"equilibrium": {"ConstantTimeLag": {"k2_delta_t": 0.025, "spin": "Free"}}}"#;
    let tides: Tides = serde_json::from_str(json).unwrap();
    assert_eq!(
        EquilibriumTide::ConstantTimeLag(test_planetary_constant_time_lag()),
        tides.equilibrium
    );
    // The spin is pseudo-synchronous unless stated.
    let json = r#"{"equilibrium": {"ConstantTimeLag": {"k2_delta_t": 0.025}}}"#;
    let tides: Tides = serde_json::from_str(json).unwrap();
    assert_eq!(
        PlanetarySpin::PseudoSynchronous,
        tides.constant_time_lag_get().unwrap().spin
    );
}

#[test]
// What is written is read back, in a self-describing format.
fn _serde_roundtrip() {
    let tides = test_stellar_tides();
    let json = serde_json::to_string(&tides).unwrap();
    assert_eq!(tides, serde_json::from_str::<Tides>(&json).unwrap());
    let json = serde_json::to_string(&Tides::default()).unwrap();
    assert_eq!(
        Tides::default(),
        serde_json::from_str::<Tides>(&json).unwrap()
    );
}

#[test]
// The keys of the former one-model-per-body layout are refused with a pointer to the new one.
fn _deserialise_former_layout_is_refused() {
    for json in [
        r#"{"StellarConstantTimeLag": {"equilibrium": {"SigmaBarStar": 1e-06}, "inertial": "FrequencyAveraged"}}"#,
        r#"{"PlanetaryConstantTimeLag": {"k2_delta_t": 0.025}}"#,
        r#"{"KaulaTides": {}}"#,
        r#""FrequencyAveraged""#,
    ] {
        let error = serde_json::from_str::<Tides>(json).unwrap_err().to_string();
        assert!(error.contains("docs/README.md"), "{error}");
    }
    // A misspelt input of the constant time lag tide is refused too.
    let json = r#"{"equilibrium": {"ConstantTimeLag": {"k2_deltat": 0.025}}}"#;
    assert!(serde_json::from_str::<Tides>(json).is_err());
}

// Checks at initialisation

#[test]
fn _initialise_star_constant_time_lag_inputs() {
    let star = |ctl: ConstantTimeLag| {
        test_particle(
            ParticleType::Star(test_star()),
            Tides {
                equilibrium: EquilibriumTide::ConstantTimeLag(ctl),
                dynamical: DynamicalTide::FrequencyAveraged,
            },
        )
    };
    // Neither input, both inputs, or a negative one.
    for ctl in [
        ConstantTimeLag::default(),
        ConstantTimeLag {
            k2_delta_t: 1e-3,
            sigma_bar_star: 1e-6,
            ..Default::default()
        },
        ConstantTimeLag {
            sigma_bar_star: -1e-6,
            ..Default::default()
        },
        ConstantTimeLag {
            k2_delta_t: f64::NAN,
            ..Default::default()
        },
    ] {
        let error = initialise_error(star(ctl));
        assert!(error.contains("sigma_bar_star"), "{error}");
    }
    // Either input alone is accepted.
    for ctl in [
        test_stellar_constant_time_lag(),
        ConstantTimeLag {
            k2_delta_t: 1e-3,
            ..Default::default()
        },
    ] {
        star(ctl).initialise(0.0).unwrap();
    }
}

#[test]
fn _initialise_planet_tides() {
    let planet = |tides: Tides| test_particle(ParticleType::Planet(test_planet()), tides);
    let with_gyration = |tides: Tides| {
        let mut planet = test_planet();
        planet.radius_of_gyration_2 = 0.254;
        test_particle(ParticleType::Planet(planet), tides)
    };
    let equilibrium = |ctl: ConstantTimeLag| Tides {
        equilibrium: EquilibriumTide::ConstantTimeLag(ctl),
        dynamical: DynamicalTide::Disabled,
    };

    // Stellar tides are refused on the planet.
    let error = initialise_error(planet(Tides {
        equilibrium: EquilibriumTide::Zahn(Zahn::default()),
        dynamical: DynamicalTide::Disabled,
    }));
    assert!(error.contains("`Zahn`"), "{error}");
    let error = initialise_error(planet(Tides {
        equilibrium: EquilibriumTide::Disabled,
        dynamical: DynamicalTide::FrequencyAveraged,
    }));
    assert!(error.contains("`FrequencyAveraged`"), "{error}");
    let error = initialise_error(planet(equilibrium(test_stellar_constant_time_lag())));
    assert!(error.contains("`sigma_bar_star`"), "{error}");

    // The time lag must be given.
    let error = initialise_error(planet(equilibrium(ConstantTimeLag::default())));
    assert!(error.contains("`k2_delta_t`"), "{error}");

    // A free spin needs the moment of inertia of the planet.
    let error = initialise_error(planet(equilibrium(test_planetary_constant_time_lag())));
    assert!(error.contains("`radius_of_gyration_2`"), "{error}");
    with_gyration(equilibrium(test_planetary_constant_time_lag()))
        .initialise(0.0)
        .unwrap();

    // A pseudo-synchronous spin does not.
    planet(equilibrium(ConstantTimeLag {
        spin: PlanetarySpin::PseudoSynchronous,
        ..test_planetary_constant_time_lag()
    }))
    .initialise(0.0)
    .unwrap();
}
