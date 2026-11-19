"""
Convert spiroid config files written before the tides were split by mechanism.

Before: one tidal model per body,
    "tides": {"StellarConstantTimeLag": {"equilibrium": {"SigmaBarStar": 1e-06},
                                         "inertial": "FrequencyAveraged"}}
    "tides": {"PlanetaryConstantTimeLag": {"k2_delta_t": 0.025, "spin": "PseudoSynchronous"}}
    "tides": {"KaulaTides": {...}}
Now: two entries per body,
    "tides": {"equilibrium": {"ConstantTimeLag": {"sigma_bar_star": 1e-06}},
              "dynamical": "FrequencyAveraged"}
    "tides": {"equilibrium": {"ConstantTimeLag": {"k2_delta_t": 0.025, "spin": "PseudoSynchronous"}},
              "dynamical": "Disabled"}
    "tides": {"equilibrium": "Disabled", "dynamical": {"KaulaTides": {...}}}

Files already in the new layout (or with tides "Disabled") are left untouched. Every other
entry of the config keeps its value; the file is rewritten with an indentation of 4.

usage: python3 convert_tides_config.py path/to/file.conf [more files or folders ...]
       python3 convert_tides_config.py --dry-run path/to/folder
A folder is searched recursively for *.conf files.
"""

import json
import sys
from pathlib import Path


def convert_tides(tides):
    """Return the tides of one body in the new layout, or None if nothing is to be changed."""
    if isinstance(tides, str):
        if tides == "Disabled":
            return None
        raise ValueError(f"unknown tides keyword {tides!r}")
    if set(tides) <= {"equilibrium", "dynamical"}:
        return None
    if len(tides) != 1:
        raise ValueError(f"unknown tides layout with keys {sorted(tides)}")

    ((model, inputs),) = tides.items()
    equilibrium = "Disabled"
    dynamical = "Disabled"
    if model == "StellarConstantTimeLag":
        old_equilibrium = inputs.get("equilibrium", "Disabled")
        if isinstance(old_equilibrium, dict) and set(old_equilibrium) == {"SigmaBarStar"}:
            equilibrium = {"ConstantTimeLag": {"sigma_bar_star": old_equilibrium["SigmaBarStar"]}}
        elif isinstance(old_equilibrium, dict) and set(old_equilibrium) == {"Zahn"}:
            equilibrium = {"Zahn": old_equilibrium["Zahn"]}
        elif old_equilibrium != "Disabled":
            raise ValueError(f"unknown stellar equilibrium tide {old_equilibrium!r}")
        old_dynamical = inputs.get("inertial", "Disabled")
        if old_dynamical not in ("Disabled", "FrequencyAveraged"):
            raise ValueError(f"unknown stellar inertial tide {old_dynamical!r}")
        dynamical = old_dynamical
    elif model == "PlanetaryConstantTimeLag":
        equilibrium = {"ConstantTimeLag": inputs}
    elif model == "KaulaTides":
        dynamical = {"KaulaTides": inputs}
    else:
        raise ValueError(f"unknown tidal model {model!r}")

    return {"equilibrium": equilibrium, "dynamical": dynamical}


def convert_file(path, dry_run):
    """Convert one config file in place. Returns True if the file needed a change."""
    config = json.loads(path.read_text())
    # A config file holds the universe under "universe"; an output snapshot is the universe.
    universe = config.get("universe", config)
    changed = False
    for body in ("central_body", "orbiting_body"):
        if body not in universe or "tides" not in universe[body]:
            continue
        converted = convert_tides(universe[body]["tides"])
        if converted is not None:
            universe[body]["tides"] = converted
            changed = True
    if changed and not dry_run:
        path.write_text(json.dumps(config, indent=4) + "\n")
    return changed


def main(arguments):
    dry_run = "--dry-run" in arguments
    targets = [Path(argument) for argument in arguments if argument != "--dry-run"]
    if not targets:
        sys.exit(__doc__)

    files = []
    for target in targets:
        files.extend(sorted(target.rglob("*.conf")) if target.is_dir() else [target])

    converted = 0
    for path in files:
        try:
            if convert_file(path, dry_run):
                converted += 1
                print(f"{'would convert' if dry_run else 'converted'}: {path}")
        except (ValueError, json.JSONDecodeError) as error:
            print(f"skipped {path}: {error}", file=sys.stderr)
    print(f"{converted} of {len(files)} files {'to convert' if dry_run else 'converted'}")


if __name__ == "__main__":
    main(sys.argv[1:])
