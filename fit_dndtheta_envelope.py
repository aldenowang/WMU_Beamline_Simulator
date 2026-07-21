#!/usr/bin/env python3
"""Fit an upper-envelope curve C * sin(theta) / sin^4(theta/2) to a
scattering-angle histogram and report the constant C.

Physics: Rutherford's differential cross-section is dsigma/dOmega ~
1/sin^4(theta/2). Converting to counts per unit *angle* (not solid angle)
requires the Jacobian dOmega = 2*pi*sin(theta)*dtheta, so

    dN/dtheta = (dN/dOmega) * dOmega/dtheta  ~  sin(theta) / sin^4(theta/2)

That extra sin(theta) factor (folded into the constant C along with the
2*pi and the beam/foil normalization) is what turns the raw 1/sin^4(theta/2)
cross-section into the count-vs-angle curve you'd overlay on a histogram
binned in theta. (analyze_scattering.py in this repo fits the bare
1/sin^4(theta/2) shape without this factor -- fine for the log-log tail
where it only rescales a straight line, but it is not the same curve as
dN/dtheta.)

Fitting method -- upper envelope, not least squares:
For each histogram bin (center theta_i, count n_i) the constant that makes
the curve pass exactly through that bin is

    C_i = n_i / f(theta_i),   f(theta) = sin(theta) / sin^4(theta/2)

Choosing C as the q-th quantile of {C_i} guarantees the curve C*f(theta)
sits at or above exactly a fraction q of the bins (n_i <= C*f(theta_i) for
every bin with C_i <= C). That's a direct, closed-form way to get a curve
under which a chosen percentage of the plotted points sit, with no
iterative optimizer needed.

Usage:
    python fit_dndtheta_envelope.py [path/to/scattering_angles_foil_*.csv]
                                     [--binwidth 1.0] [--quantile 0.95]
                                     [--plot]

With no CSV argument, uses the newest scattering_angles_foil_*.csv in data/.
Prints C and writes a <csv-stem>_envelope.json with the fitted curve
(for overlaying on other plots, e.g. the HTML histogram artifact) and,
with --plot, a <csv-stem>_envelope.png.
"""

import argparse
import csv
import glob
import json
import math
import os
import sys

import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

COLOR_BARS = "#2a78d6"
COLOR_CURVE = "#eb6834"      # 95th-percentile envelope
COLOR_CURVE_100 = "#1baf7a"  # 100th-percentile (true upper bound) envelope
COLOR_MISS = "#e34948"
COLOR_INK = "#0b0b0b"
COLOR_SECONDARY_INK = "#52514e"
COLOR_GRID = "#e1e0d9"
COLOR_AXIS = "#c3c2b7"


def find_latest_foil_csv():
    candidates = sorted(glob.glob(os.path.join(DATA_DIR, "scattering_angles_foil_*.csv")))
    if not candidates:
        raise SystemExit(f"No scattering_angles_foil_*.csv files found in {DATA_DIR}")
    return candidates[-1]


def load_angles(csv_path):
    angles = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            angles.append(float(row["angle_deg"]))
    return np.array(angles)


def dndtheta_shape(theta_deg):
    """f(theta) = sin(theta) / sin^4(theta/2), theta in degrees. Diverges as
    theta -> 0, which is expected: Rutherford single-scattering only governs
    the sparse large-angle tail -- the near-zero bulk is multiple scattering
    and will naturally sit far *under* this curve, not on it."""
    theta_rad = np.radians(theta_deg)
    return np.sin(theta_rad) / np.sin(theta_rad / 2.0) ** 4


def fit_envelope(angles, binwidth=1.0, quantile=0.95):
    max_angle = angles.max()
    nbins = int(math.ceil(max_angle / binwidth))
    edges = np.arange(nbins + 1) * binwidth
    counts, edges = np.histogram(angles, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    mask = counts > 0
    f_vals = dndtheta_shape(centers[mask])
    c_i = counts[mask] / f_vals
    C = float(np.quantile(c_i, quantile))

    predicted = C * dndtheta_shape(centers[mask])
    coverage = float(np.mean(counts[mask] <= predicted))

    return {
        "C": C,
        "quantile_requested": quantile,
        "coverage_achieved": coverage,
        "n_bins_fit": int(mask.sum()),
        "binwidth": binwidth,
        "centers": centers,
        "counts": counts,
        "mask": mask,
    }


def style_axes(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=COLOR_GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(COLOR_AXIS)
    ax.tick_params(colors=COLOR_SECONDARY_INK)


def plot_envelope(result, angles, out_path):
    import matplotlib.pyplot as plt

    centers, counts, mask = result["centers"], result["counts"], result["mask"]
    C = result["C"]

    theta_grid = np.linspace(max(centers[mask].min() * 0.3, 1e-2), centers[mask].max(), 400)
    curve = C * dndtheta_shape(theta_grid)

    above = mask & (counts > C * dndtheta_shape(centers))

    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    style_axes(ax)

    ax.bar(centers[mask & ~above], counts[mask & ~above], width=result["binwidth"] * 0.92,
           color=COLOR_BARS, edgecolor="#fcfcfb", linewidth=0.3,
           label=f"Histogram bins under envelope (n={angles.size})")
    if above.any():
        ax.bar(centers[above], counts[above], width=result["binwidth"] * 0.92,
               color=COLOR_MISS, edgecolor="#fcfcfb", linewidth=0.3,
               label=f"Bins above envelope ({above.sum()}, {100*(1-result['coverage_achieved']):.1f}%)")

    ax.plot(theta_grid, curve, color=COLOR_CURVE, linewidth=2,
            label=rf"Envelope: $dN/d\theta = {C:.3g} \cdot \sin\theta / \sin^4(\theta/2)$"
                  f"\n({result['quantile_requested']*100:.0f}th pct fit, "
                  f"{result['coverage_achieved']*100:.1f}% of bins covered)")

    ax.set_yscale("log")
    ax.set_xlabel("Scattering angle from beam path, θ (degrees)", color=COLOR_INK)
    ax.set_ylabel("Number of particles per bin (log scale)", color=COLOR_INK)
    ax.set_title("dN/dθ upper-envelope fit to foil-exit scattering angles", color=COLOR_INK, fontsize=13)
    ax.legend(frameon=False, labelcolor=COLOR_INK, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    return out_path


def plot_dual_envelope(result_95, result_100, angles, out_path):
    import matplotlib.pyplot as plt

    centers, counts, mask = result_100["centers"], result_100["counts"], result_100["mask"]

    theta_grid = np.linspace(max(centers[mask].min() * 0.3, 1e-2), centers[mask].max(), 400)
    curve_95 = result_95["C"] * dndtheta_shape(theta_grid)
    curve_100 = result_100["C"] * dndtheta_shape(theta_grid)

    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    style_axes(ax)

    ax.bar(centers[mask], counts[mask], width=result_100["binwidth"] * 0.92,
           color=COLOR_BARS, edgecolor="#fcfcfb", linewidth=0.3,
           label=f"Simulated foil-exit angles (n={angles.size})")

    ax.plot(theta_grid, curve_95, color=COLOR_CURVE, linewidth=2, linestyle="-",
            label=rf"95th-pct envelope: $C={result_95['C']:.3g}$"
                  f" ({result_95['coverage_achieved']*100:.1f}% of bins covered)")
    ax.plot(theta_grid, curve_100, color=COLOR_CURVE_100, linewidth=2, linestyle="--",
            label=rf"100th-pct envelope (true upper bound): $C={result_100['C']:.3g}$")

    ax.set_yscale("log")
    ax.set_xlabel("Scattering angle from beam path, θ (degrees)", color=COLOR_INK)
    ax.set_ylabel("Number of particles per bin (log scale)", color=COLOR_INK)
    ax.set_title("dN/dθ upper-envelope fits to foil-exit scattering angles", color=COLOR_INK, fontsize=13)
    ax.legend(frameon=False, labelcolor=COLOR_INK, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="?", help="scattering_angles_foil_*.csv (default: newest in data/)")
    parser.add_argument("--binwidth", type=float, default=1.0, help="histogram bin width in degrees (default 1.0)")
    parser.add_argument("--quantile", type=float, default=0.95, help="fraction of bins to sit under the curve (default 0.95)")
    parser.add_argument("--plot", action="store_true", help="also write a PNG plot")
    args = parser.parse_args()

    csv_path = args.csv or find_latest_foil_csv()
    angles = load_angles(csv_path)
    if angles.size == 0:
        raise SystemExit(f"{csv_path} has no rows")

    result = fit_envelope(angles, binwidth=args.binwidth, quantile=args.quantile)

    print(f"Read {angles.size} angles from {csv_path}")
    print(f"dN/dtheta = C * sin(theta) / sin^4(theta/2)")
    print(f"  C                 = {result['C']:.6g}")
    print(f"  quantile fit      = {result['quantile_requested']*100:.1f}%")
    print(f"  coverage achieved = {result['coverage_achieved']*100:.2f}% of {result['n_bins_fit']} non-empty bins")

    q_tag = f"q{args.quantile*100:g}".replace(".", "p")
    out_json = os.path.splitext(csv_path)[0] + f"_envelope_{q_tag}.json"
    theta_grid = np.linspace(max(result["centers"][result["mask"]].min() * 0.3, 1e-2),
                              result["centers"][result["mask"]].max(), 200)
    curve = result["C"] * dndtheta_shape(theta_grid)
    with open(out_json, "w") as f:
        json.dump({
            "formula": "dN/dtheta = C * sin(theta_deg) / sin(theta_deg/2)^4, theta in degrees",
            "C": result["C"],
            "quantile_requested": result["quantile_requested"],
            "coverage_achieved": result["coverage_achieved"],
            "binwidth": result["binwidth"],
            "curve_theta_deg": theta_grid.tolist(),
            "curve_dndtheta": curve.tolist(),
        }, f, indent=2)
    print(f"Wrote curve data to {out_json}")

    if args.plot:
        result_95 = result if args.quantile == 0.95 else fit_envelope(angles, binwidth=args.binwidth, quantile=0.95)
        result_100 = result if args.quantile == 1.0 else fit_envelope(angles, binwidth=args.binwidth, quantile=1.0)
        out_png = os.path.splitext(csv_path)[0] + "_envelope.png"
        plot_dual_envelope(result_95, result_100, angles, out_png)
        print(f"Wrote plot to {out_png}")


if __name__ == "__main__":
    main()
