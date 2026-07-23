#!/usr/bin/env python3
"""Plot simulated alpha particle scattering angles against the Rutherford 1/sin^4(theta/2) prediction.

Reads the scattering_angles_foil_*.csv files written by RunAction (one
"angle_deg" column, one row per primary particle): the angle the instant
each primary leaves the silicon foil -- the actual scattering event, spanning
the full angular range including the rare large-angle single-Coulomb-
scattering tail Rutherford's law predicts.

Usage:
    python analyze_scattering.py [path/to/scattering_angles_foil_*.csv]

With no argument, uses the newest foil CSV in data/.
"""

import csv
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# dataviz reference palette: categorical slot 1 (blue, bars) and slot 6
# (orange, theory curve) -- validated as a CVD-safe pair.
COLOR_BARS = "#2a78d6"
COLOR_CURVE = "#eb6834"
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


def rutherford_shape(theta_deg):
    """Unnormalized Rutherford differential cross-section shape, 1/sin^4(theta/2)."""
    theta_rad = np.radians(theta_deg)
    return 1.0 / np.sin(theta_rad / 2.0) ** 4


def bin_integrated_theory(edges, n_sub=50):
    """Expected relative counts per bin: integral of rutherford_shape over
    each [edges[i], edges[i+1]], via fine-grained trapezoidal quadrature
    (bins can be log-wide, so a midpoint estimate isn't accurate enough)."""
    totals = np.empty(len(edges) - 1)
    for i in range(len(edges) - 1):
        x = np.linspace(edges[i], edges[i + 1], n_sub)
        totals[i] = np.trapezoid(rutherford_shape(x), x)
    return totals


def style_axes(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=COLOR_GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(COLOR_AXIS)
    ax.tick_params(colors=COLOR_SECONDARY_INK)


def plot_foil(angles, out_path):
    """Full-range foil-exit angles: log-spaced bins (the data spans several
    orders of magnitude in angle) on a log-log plot, the classic way to see
    Rutherford's power-law tail as a straight line.

    The near-zero-to-few-degree bulk is dominated by *multiple* small-angle
    Coulomb scattering (many tiny deflections accumulating to a roughly
    Gaussian spread), not single scattering off one nucleus -- that's why the
    histogram rises to a peak and then falls, instead of falling
    monotonically like 1/sin^4(theta/2). Rutherford's law only governs the
    sparse large-angle tail beyond that peak, where a single hard scattering
    event dominates over the cumulative multiple-scattering background --
    exactly why the historical Geiger-Marsden experiment only ever counted
    flashes at large angles and never tried to match the small-angle
    behaviour. So the fit here is deliberately restricted to bins at/after
    the observed peak, and the theory curve is only drawn over that same
    domain -- fitting or drawing it over the multiple-scattering-dominated
    region would be physically meaningless.
    """
    theta_min = max(angles.min(), 1e-3)
    edges = np.logspace(np.log10(theta_min), np.log10(angles.max()), 45)
    counts, edges = np.histogram(angles, bins=edges)
    centers = np.sqrt(edges[:-1] * edges[1:])  # geometric center for log bins
    widths = np.diff(edges)

    theory_bins = bin_integrated_theory(edges)

    peak_idx = int(np.argmax(counts))
    tail = np.zeros(len(counts), dtype=bool)
    tail[peak_idx:] = True
    tail &= counts > 0

    scale = np.sum(counts[tail] * theory_bins[tail]) / np.sum(theory_bins[tail] ** 2)
    theory_counts = scale * theory_bins

    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    style_axes(ax)

    ax.bar(
        centers,
        counts,
        width=widths * 0.92,
        color=COLOR_BARS,
        edgecolor="#fcfcfb",
        linewidth=0.3,
        label=f"Simulated foil-exit angles (n={angles.size})",
    )
    ax.plot(
        centers[peak_idx:],
        theory_counts[peak_idx:],
        color=COLOR_CURVE,
        linewidth=2,
        marker="o",
        markersize=3,
        label=r"Rutherford prediction $\propto 1/\sin^4(\theta/2)$"
        "\n(fit to large-angle tail past the multiple-scattering peak)",
    )
    ax.axvline(centers[peak_idx], color=COLOR_AXIS, linewidth=1, linestyle="--")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Scattering angle from beam path, θ (degrees, log scale)", color=COLOR_INK)
    ax.set_ylabel("Number of particles (log scale)", color=COLOR_INK)
    ax.set_title(
        "Rutherford scattering: alpha particle angle on exiting the 100 µm silicon foil",
        color=COLOR_INK,
        fontsize=13,
    )
    ax.legend(frameon=False, labelcolor=COLOR_INK, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    return out_path


def main():
    foil_csv = sys.argv[1] if len(sys.argv) > 1 else find_latest_foil_csv()

    foil_angles = load_angles(foil_csv)
    if foil_angles.size == 0:
        raise SystemExit(f"{foil_csv} has no rows")
    foil_out = plot_foil(foil_angles, os.path.splitext(foil_csv)[0] + "_plot.png")
    print(f"Read {foil_angles.size} angles from {foil_csv}")
    print(f"Wrote plot to {foil_out}")


if __name__ == "__main__":
    main()
