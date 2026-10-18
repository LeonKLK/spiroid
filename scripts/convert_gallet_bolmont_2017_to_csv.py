"""
Builds a one-zone (no radiative core) spiroid star file from a Gallet & Bolmont 2017 track,
reproducing the rigid star of Posidonius' `GalletBolmont2017` evolution: the radius follows
the track and the moment of inertia is I = rg2 M R(t)^2 with a constant rg2 (0.07 by default).

The radiative zone is emptied (radiative_mass = radiative_moment_of_inertia = 0,
convective_radius = 0, i.e. no core), so the core-envelope coupling and the mass transfer
torques vanish and the envelope equation reduces to dJ/dt = T_wind + T_tide.

Input, one of:
    --track M_10_Z_0134.dat   Posidonius input/Gallet_Bolmont_2017 file (space separated,
                              optional header): column 1 log10(age / yr), column 3 T_eff (K),
                              column 4 R / R_sun (the columns Posidonius reads).
    --history history.txt     Posidonius history (tab separated, diagnostic.py columns). The
                              radius Posidonius used is piecewise linear in time between the
                              track nodes; the nodes are recovered from the slope breaks of the
                              `stellar_radius` column and written out, so the output reproduces
                              the radius of that run exactly (to ~1e-10 R_sun).
                              `--initial-age` (yr) is the track age of the history's time 0.

Output columns (StarCsv, src/universe/particles/star/star_csv.rs): age (yr), radius (R_sun),
mass (M_sun), convective_radius (R_sun), radiative_mass (M_sun), radiative_moment_of_inertia
and convective_moment_of_inertia (M R^2), luminosity (L_sun).

`--mass-kg` sets the stellar mass written to the file (default 1.9818e30 kg, the Posidonius
M_SUN) so that G M, and hence the mean motion, matches Posidonius.

Examples:
    python3 scripts/convert_gallet_bolmont_2017_to_csv.py --track input/Gallet_Bolmont_2017/M_10_Z_0134.dat \
        examples/data/star/evolution/gallet_bolmont_2017_M_10_Z_0134_adapted.csv
    python3 scripts/convert_gallet_bolmont_2017_to_csv.py --history sma_03_kaula_outward_history.txt \
        --initial-age 5e6 examples/data/star/evolution/gallet_bolmont_2017_M_10_Z_0134_adapted.csv
"""

import argparse
import sys

import numpy as np

sys.dont_write_bytecode = True

# spiroid's (astro-const) solar mass, the unit of the `mass` column.
SOLAR_MASS = 1.327_124_420_99e20 / 6.674_28e-11
POSIDONIUS_SOLAR_MASS = 1.9818e30
SOLAR_EFFECTIVE_TEMPERATURE = 5772.0


def read_track(path):
    ages, radii, teffs = [], [], []
    with open(path) as f:
        for line in f:
            fields = line.split()
            try:
                log_age, teff, radius = float(fields[0]), float(fields[2]), float(fields[3])
            except (ValueError, IndexError):
                continue  # header
            ages.append(10.0**log_age)
            teffs.append(teff)
            radii.append(radius)
    radii = np.array(radii)
    luminosity = radii**2 * (np.array(teffs) / SOLAR_EFFECTIVE_TEMPERATURE) ** 4
    return np.array(ages), radii, luminosity


def read_history(path, initial_age):
    with open(path) as f:
        header = f.readline().split()
    columns = (header.index("current_time"), header.index("stellar_radius"))
    t, radius = np.loadtxt(path, skiprows=1, usecols=columns, unpack=True)
    # Slope breaks: the radius is linear between track nodes. A node falls inside one sampled
    # interval, whose slope then differs from both neighbours.
    slope = np.diff(radius) / np.diff(t)
    breaks = np.where(np.abs(np.diff(slope)) * np.median(np.diff(t)) > 1e-11)[0]
    groups = []
    for j in breaks:
        if groups and j - groups[-1][-1] <= 1:
            groups[-1].append(j)
        else:
            groups.append([j])
    # Clean (single line) stretches of samples between the breaks.
    segments, start = [], 0
    for group in groups:
        segments.append((start, group[0] + 1))
        start = group[-1] + 2
    segments.append((start, len(t) - 1))
    lines = []
    for a, b in segments:
        b = max(b, a + 1)
        m = (radius[b] - radius[a]) / (t[b] - t[a])
        lines.append((m, radius[a] - m * t[a]))
    # Nodes at the intersections of consecutive lines.
    node_t, node_r = [t[0]], [radius[0]]
    for (m1, c1), (m2, c2) in zip(lines[:-1], lines[1:]):
        x = (c2 - c1) / (m1 - m2)
        node_t.append(x)
        node_r.append(m1 * x + c1)
    node_t.append(t[-1])
    node_r.append(radius[-1])
    node_t, node_r = np.array(node_t), np.array(node_r)
    error = np.abs(np.interp(t, node_t, node_r) - radius).max()
    print(f"{len(node_t)} track nodes recovered, max radius error {error:.2e} R_sun")
    # Luminosity is not used without magnetism; any positive value will do.
    return node_t + initial_age, node_r, np.ones_like(node_r)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--track", help="Gallet & Bolmont 2017 .dat file")
    source.add_argument("--history", help="Posidonius history file")
    parser.add_argument("--initial-age", type=float, default=5e6, help="track age (yr) of the history's time 0")
    parser.add_argument("--rg2", type=float, default=0.07, help="constant squared radius of gyration")
    parser.add_argument("--mass-kg", type=float, default=POSIDONIUS_SOLAR_MASS, help="stellar mass (kg)")
    parser.add_argument("output", help="spiroid star CSV")
    args = parser.parse_args()

    if args.track:
        ages, radii, luminosity = read_track(args.track)
    else:
        ages, radii, luminosity = read_history(args.history, args.initial_age)

    zeros = np.zeros_like(radii)
    table = np.column_stack(
        [
            ages,
            radii,
            zeros + args.mass_kg / SOLAR_MASS,
            zeros,  # convective_radius: base of the envelope, no core
            zeros,  # radiative_mass
            zeros,  # radiative_moment_of_inertia
            zeros + args.rg2,  # convective_moment_of_inertia (M R^2)
            luminosity,
        ]
    )
    header = "age,radius,mass,convective_radius,radiative_mass,radiative_moment_of_inertia,convective_moment_of_inertia,luminosity"
    np.savetxt(args.output, table, delimiter=",", header=header, comments="", fmt="%.18e")
    print(f"wrote {len(ages)} rows to {args.output}")


if __name__ == "__main__":
    main()
