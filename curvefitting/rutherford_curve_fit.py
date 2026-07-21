#!/usr/bin/env python3
"""Least-squares fit of the Rutherford dN/dtheta curve to a simulated
scattering-angle histogram, with the theoretical amplitude computed
independently from first principles so the two can be compared directly.

Physics recap
--------------
Rutherford's differential cross-section (point-charge, non-relativistic):

    dsigma/dOmega = ( Z1 * Z2 * e^2 / (4*E) )^2 * 1 / sin^4(theta/2)      (1)

e^2 here means e^2/(4*pi*eps0) = 1.44 MeV*fm, E is the incident kinetic
energy, theta is the CM/lab scattering angle. This is a cross-section per
target nucleus per steradian -- to get *counts* you still need to multiply
by how many incident particles there were and how many target nuclei per
unit area they could have scattered off:

    dN/dOmega = N_total * n_areal * dsigma/dOmega                        (2)

    n_areal = rho * N_A * t / M     (target nuclei / cm^2, thin-foil)

Histograms bin in theta, not Omega, so converting (2) to a per-angle count
needs the solid-angle Jacobian dOmega = 2*pi*sin(theta) dtheta (theta in
RADIANS -- this is where the sin(theta)/sin^4(theta/2) shape comes from):

    dN/dtheta = 2*pi * N_total * n_areal * (Z1*Z2*e^2/4E)^2
                * sin(theta) / sin^4(theta/2)
              = C_theory * sin(theta) / sin^4(theta/2)                   (3)

Two unit traps that silently break the amplitude match between a fitted C
and this C_theory:

  1. sin() in (1)/(3) must see RADIANS. If the histogram's theta axis is
     in degrees and degree-values get passed straight into np.sin(), the
     shape (and therefore the fitted amplitude) is wrong -- np.sin(5.0)
     treats 5.0 as radians, not 5 degrees.
  2. dN/dtheta is a *density* (counts per unit angle), not a raw bin
     count. A raw histogram count in a bin of width Delta-theta estimates
     dN/dtheta * Delta-theta, not dN/dtheta itself. If Delta-theta isn't
     exactly 1 (in whatever angle unit the fit function's theta argument
     is in), the fitted amplitude is off by that bin-width factor. Since
     the theoretical curve (3) is written in radians, the histogram must
     be normalized as counts / Delta-theta_radians to be compared against
     it on equal footing -- normalizing by the degree bin width instead
     leaves a stray (pi/180) factor in the fitted amplitude.

This script converts to radians and normalizes by the radian bin width
before fitting, so C_fit should land close to C_theory once those two
issues are controlled for. Any *remaining* gap is real physics/geometry
(nuclear size, nuclear or electron screening, nonzero detector angular
resolution, foil straggling, etc.) rather than a units bug.

Usage:
    python rutherford_curve_fit.py [path/to/scattering_angles_foil_*.csv]
                                    [--binwidth 1.0] [--theta-min deg]
                                    [--Z1 1] [--Z2 79] [--E 6.0]
                                    [--rho 19.3] [--t 0.0001] [--M 197.0]

With no CSV argument, uses the newest scattering_angles_foil_*.csv in
../data relative to this file.
"""

import argparse
import csv
import glob
import os

import numpy as np
from scipy.optimize import curve_fit

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")

# dataviz reference palette (validated CVD-safe triple, see
# fit_dndtheta_envelope.py / MEMORY for the validation run):
COLOR_BARS = "#2a78d6"      # histogram
COLOR_THEORY = "#eb6834"    # theoretical curve (computed, no fitting)
COLOR_FIT = "#1baf7a"       # best-fit curve
COLOR_INK = "#0b0b0b"
COLOR_SECONDARY_INK = "#52514e"
COLOR_GRID = "#e1e0d9"
COLOR_AXIS = "#c3c2b7"

E2_MEV_FM = 1.44      # e^2 / (4 pi eps0), MeV*fm  (Z1=Z2=1 point charges)
FM2_TO_CM2 = 1e-26
N_A = 6.022e23        # Avogadro's number, 1/mol


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


# ---------------------------------------------------------------------
# Model function used by curve_fit. theta comes in DEGREES (matching the
# histogram's native x-axis) and is converted to radians *inside* the
# model, before the only two sin() calls -- this is fix #1 from the
# docstring above.
# ---------------------------------------------------------------------
def dndtheta_model(theta_deg, C):
    theta_rad = np.radians(theta_deg)
    return C * np.sin(theta_rad) / np.sin(theta_rad / 2.0) ** 4


def theoretical_amplitude(n_total, Z1, Z2, E_mev, rho, t_cm, M):
    """C_theory such that dN/dtheta = C_theory * sin(theta)/sin^4(theta/2),
    theta in radians. See docstring eq. (3)."""
    C1_fm2 = (Z1 * Z2 * E2_MEV_FM / (4.0 * E_mev)) ** 2
    C1_cm2 = C1_fm2 * FM2_TO_CM2
    n_areal = rho * N_A * t_cm / M  # target nuclei / cm^2
    C2 = n_total * n_areal
    C_theory = 2 * np.pi * C1_cm2 * C2
    return C_theory, {"C1_cm2": C1_cm2, "n_areal_per_cm2": n_areal, "C2": C2}


def style_axes(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=COLOR_GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(COLOR_AXIS)
    ax.tick_params(colors=COLOR_SECONDARY_INK)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="?", help="scattering_angles_foil_*.csv (default: newest in data/)")
    parser.add_argument("--binwidth", type=float, default=1.0, help="histogram bin width in degrees (default 1.0)")
    parser.add_argument("--theta-min", type=float, default=None,
                         help="degrees; restrict the fit to theta >= this (default: auto-detect the "
                              "histogram peak and fit only the tail past it, since Rutherford's "
                              "point-charge formula doesn't describe the multiple-scattering bulk)")
    parser.add_argument("--Z1", type=float, default=1.0, help="projectile charge number (default 1, proton)")
    parser.add_argument("--Z2", type=float, default=79.0, help="target charge number (default 79, gold)")
    parser.add_argument("--E", type=float, default=6.0, help="incident kinetic energy, MeV (default 6.0)")
    parser.add_argument("--rho", type=float, default=19.3, help="target density, g/cm^3 (default 19.3, gold)")
    parser.add_argument("--t", type=float, default=0.0001, help="foil thickness, cm (default 0.0001 = 1 um, matches DetectorConstruction.cc kFoilHalfZ)")
    parser.add_argument("--M", type=float, default=197.0, help="target molar mass, g/mol (default 197, gold)")
    parser.add_argument("--plot", action="store_true", help="also write a PNG plot")
    parser.add_argument("--fit-only", action="store_true",
                         help="in the plot, show only the histogram + best-fit curve (omit the "
                              "theoretical/computed curve and its legend entry)")
    args = parser.parse_args()

    csv_path = args.csv or find_latest_foil_csv()
    angles = load_angles(csv_path)
    n_total = angles.size
    if n_total == 0:
        raise SystemExit(f"{csv_path} has no rows")

    # ---- histogram -----------------------------------------------------
    nbins = int(np.ceil(angles.max() / args.binwidth))
    edges = np.arange(nbins + 1) * args.binwidth
    counts, edges = np.histogram(angles, bins=edges)
    centers_deg = 0.5 * (edges[:-1] + edges[1:])

    # fix #2: normalize by the bin width IN RADIANS, since the model (and
    # C_theory) are both derived with dOmega = 2*pi*sin(theta) dtheta_rad.
    binwidth_rad = np.radians(args.binwidth)
    density = counts / binwidth_rad  # dN/dtheta estimate, per radian

    theta_min = args.theta_min
    if theta_min is None:
        peak_idx = int(np.argmax(counts))
        theta_min = centers_deg[peak_idx]
    mask = (counts > 0) & (centers_deg >= theta_min)
    if mask.sum() < 2:
        raise SystemExit(f"Only {mask.sum()} non-empty bins at theta >= {theta_min} deg -- can't fit")

    # ---- fit -------------------------------------------------------------
    # Poisson counting error on each bin (sqrt(N)), propagated through the
    # counts -> density normalization. Fitting *unweighted* would let the
    # handful of loudest near-peak bins dominate the residual sum and
    # starve the power-law tail of any influence on the fit -- exactly the
    # kind of thing that produces an amplitude that "doesn't match theory"
    # for reasons that have nothing to do with units.
    sigma = np.sqrt(counts[mask]) / binwidth_rad
    p0 = [density[mask].max()]
    popt, pcov = curve_fit(
        dndtheta_model, centers_deg[mask], density[mask],
        p0=p0, sigma=sigma, absolute_sigma=True, bounds=(0, np.inf),
    )
    C_fit = popt[0]
    C_fit_err = float(np.sqrt(pcov[0, 0]))

    # ---- goodness of fit -------------------------------------------------
    # n_bins = number of non-empty bins in the fit region (theta >= theta_min);
    # n_params = 1 (just the amplitude C -- theta_min itself is a fixed cut,
    # not a fitted parameter, and Z1/Z2/E/rho/t/M only enter C_theory, not
    # this fit). dof = n_bins - n_params.
    n_bins_fit = int(mask.sum())
    n_params = len(popt)
    dof = n_bins_fit - n_params
    residuals = (density[mask] - dndtheta_model(centers_deg[mask], C_fit)) / sigma
    chi2 = float(np.sum(residuals ** 2))
    reduced_chi2 = chi2 / dof if dof > 0 else float("nan")
    from scipy.stats import chi2 as chi2_dist
    p_value = float(chi2_dist.sf(chi2, dof)) if dof > 0 else float("nan")

    # ---- theory ------------------------------------------------------
    C_theory, parts = theoretical_amplitude(n_total, args.Z1, args.Z2, args.E, args.rho, args.t, args.M)

    print(f"Read {n_total} angles from {csv_path}")
    print(f"Histogram: {args.binwidth} deg bins, fit restricted to theta >= {theta_min:.3f} deg "
          f"({n_bins_fit} non-empty bins)")
    print()
    print("Theoretical amplitude  (dN/dtheta = C_theory * sin(theta)/sin^4(theta/2), theta in radians):")
    print(f"  C1 = (Z1*Z2*e^2/4E)^2       = {parts['C1_cm2']:.6e} cm^2")
    print(f"  n_areal = rho*N_A*t/M       = {parts['n_areal_per_cm2']:.6e} nuclei/cm^2")
    print(f"  C2 = N_total * n_areal      = {parts['C2']:.6e}")
    print(f"  C_theory = 2*pi * C1 * C2   = {C_theory:.6e}")
    print()
    print("Fitted amplitude (scipy.optimize.curve_fit, Poisson-weighted, theta-domain restricted to tail):")
    print(f"  C_fit = {C_fit:.6e} +/- {C_fit_err:.2e}")
    print()
    print(f"ratio C_fit / C_theory = {C_fit / C_theory:.6g}")
    print()
    print("Chi-square goodness of fit (Poisson-weighted residuals, fit region only):")
    print(f"  n_bins (fit region, theta >= {theta_min:.3f} deg) = {n_bins_fit}")
    print(f"  n_params (free parameters: C only)                = {n_params}")
    print(f"  dof = n_bins - n_params                            = {dof}")
    print(f"  chi2                                               = {chi2:.4f}")
    print(f"  reduced chi2 = chi2 / dof                          = {reduced_chi2:.4f}")
    print(f"  p-value (chi2.sf, dof={dof})                       = {p_value:.4g}")

    if args.plot:
        import matplotlib.pyplot as plt

        theta_plot = np.linspace(centers_deg[mask].min(), centers_deg.max(), 400)
        fit_curve = dndtheta_model(theta_plot, C_fit)

        fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
        fig.patch.set_facecolor("#fcfcfb")
        style_axes(ax)

        ax.bar(centers_deg[counts > 0], density[counts > 0], width=args.binwidth * 0.92,
               color=COLOR_BARS, edgecolor="#fcfcfb", linewidth=0.3,
               label=f"Simulated foil-exit angles, dN/d$\\theta$ (n={n_total})")
        if not args.fit_only:
            theory_curve = dndtheta_model(theta_plot, C_theory)
            ax.plot(theta_plot, theory_curve, color=COLOR_THEORY, linewidth=2, linestyle="--",
                    label=f"Theoretical (computed): $C_{{theory}}$ = {C_theory:.3e}")
        ax.plot(theta_plot, fit_curve, color=COLOR_FIT, linewidth=2.5, zorder=5,
                label=f"Best fit (curve_fit): $C_{{fit}}$ = {C_fit:.3e}")
        ax.axvline(theta_min, color=COLOR_AXIS, linewidth=1, linestyle=":")

        ax.set_yscale("log")
        ax.set_xlabel(r"Scattering angle, $\theta$ (degrees)", color=COLOR_INK)
        ax.set_ylabel(r"$dN/d\theta$ (counts per radian, log scale)", color=COLOR_INK)
        title = ("Rutherford dN/d$\\theta$: best fit vs. simulated histogram" if args.fit_only
                 else "Rutherford dN/d$\\theta$: theory vs. best fit vs. simulated histogram")
        ax.set_title(title, color=COLOR_INK, fontsize=13)
        ax.legend(frameon=False, labelcolor=COLOR_INK, fontsize=9)

        eqn_text = (
            r"$\dfrac{dN}{d\theta} = C \cdot \dfrac{\sin\theta}{\sin^4(\theta/2)}$"
            "\n"
            fr"$C = {C_fit:.3f} \pm {C_fit_err:.3f}$ counts/rad"
            "\n"
            fr"dof = {n_bins_fit} bins $-$ {n_params} = {dof}"
        )
        ax.text(0.60, 0.60, eqn_text, transform=ax.transAxes, fontsize=11,
                color=COLOR_INK, ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="#fcfcfb",
                          edgecolor=COLOR_FIT, linewidth=1.2))

        fig.tight_layout()

        out_png = os.path.splitext(csv_path)[0] + "_rutherford_fit.png"
        fig.savefig(out_png, facecolor=fig.get_facecolor())
        print(f"Wrote plot to {out_png}")


if __name__ == "__main__":
    main()
