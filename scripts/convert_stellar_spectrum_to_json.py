"""
Converts a tabulated love number spectrum of the stellar dynamical tide (inertial waves in
the convective envelope) into the `SpectrumFile` JSON read by spiroid for the
`StarConvectiveEnvelope` kaula composition (see src/universe/effects/tides/kaula/love_number.rs).

Input: a whitespace separated text file with one row per tidal frequency and three columns,
    omega        tidal frequency (rad.s-1) at the rotation period the spectrum was computed for,
    Im(k2)       imaginary part of the l=m=2 love number (sign of omega),
    Re(k2)       real part of the l=m=2 love number (hydrostatic, constant).

Output: {"spin_spec": Omega (rad.s-1), "spectrum": {"Interpolate1D": {"info", "x_vals", "data"}}}
    x_vals are the dimensionless tidal frequencies omega / Omega (spiroid multiplies them back
    by spin_spec when the file is loaded), data is [Re(k2), Im(k2)].

Example:
    python3 scripts/convert_stellar_spectrum_to_json.py \
        examples/data/star/tides/kaula/alpha0.51600_P1p2_Ek1p5em6.txt \
        examples/data/star/tides/kaula/spectrum_alpha0.516_P1.2d_Ek1.5e-6.json \
        --period 1.2 --info "1 Msun, alpha = 0.516, Ek = 1.5e-6, from A. Astoul"
"""

import argparse
import json
import math
import sys

sys.dont_write_bytecode = True

from units import SECONDS_IN_DAY


def read_spectrum(path):
    omegas, re_k2, im_k2 = [], [], []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) != 3:
                continue
            omega, im, re = (float(x) for x in fields)
            omegas.append(omega)
            im_k2.append(im)
            re_k2.append(re)
    return omegas, re_k2, im_k2


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="text spectrum: omega (rad/s), Im(k2), Re(k2)")
    parser.add_argument("output", help="output json path")
    parser.add_argument("--period", type=float, required=True, help="rotation period of the spectrum (days)")
    parser.add_argument("--info", default="", help="free text stored in the file")
    args = parser.parse_args()

    spin_spec = 2.0 * math.pi / (args.period * SECONDS_IN_DAY)
    omegas, re_k2, im_k2 = read_spectrum(args.input)

    # The interpolator needs strictly increasing frequencies.
    rows = sorted(zip(omegas, re_k2, im_k2))
    x_vals = [omega / spin_spec for omega, _, _ in rows]
    data = [[re, im] for _, re, im in rows]
    if any(b <= a for a, b in zip(x_vals, x_vals[1:])):
        sys.exit("tidal frequencies must be strictly increasing")

    info = (
        "Complex l=m=2 love number spectrum of the stellar dynamical tide (inertial waves). "
        f"x_vals: dimensionless tidal frequency omega / spin_spec; spin_spec = {spin_spec!r} rad/s "
        f"(P = {args.period} days). data: [Re(k2), Im(k2)]. Converted from {args.input}. {args.info}"
    ).strip()

    spectrum = {
        "spin_spec": spin_spec,
        "spectrum": {"Interpolate1D": {"info": info, "x_vals": x_vals, "data": data}},
    }
    with open(args.output, "w") as f:
        json.dump(spectrum, f)
    print(f"wrote {args.output}: {len(x_vals)} frequencies in [{x_vals[0]:.4g}, {x_vals[-1]:.4g}] x spin_spec")


if __name__ == "__main__":
    main()
