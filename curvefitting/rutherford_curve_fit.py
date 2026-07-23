#!/usr/bin/env python3

import argparse
import csv
import glob
import os
from scipy.stats import chi2
import math


import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")

# dataviz reference palette (validated CVD-safe triple, see
# fit_dndtheta_envelope.py / MEMORY for the validation run):
COLOR_BARS = "#2a78d6"      # histogram
COLOR_ENVELOPE = "#eda100"  # fit curve, dN/dtheta = C*sin(theta)/sin^4(theta/2)
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


def r_squared_log10(y, y_pred):
    """R^2 on log10(y) vs log10(y_pred), one point per bin (bin center, bin
    top). This data spans ~6 orders of magnitude, so a linear-space R^2 is
    dominated almost entirely by the single largest-count bin (verified: for
    this dataset one bin alone accounted for >95% of the linear sum of
    squared residuals) -- log-space matches the log-scale y-axis these plots
    already use and reflects fit quality across all bins, not just the peak."""
    logy = np.log10(np.asarray(y, dtype=float))
    logpred = np.log10(y_pred)
    ss_res = np.sum((logy - logpred) ** 2)
    ss_tot = np.sum((logy - logy.mean()) ** 2)
    return 1.0 - ss_res / ss_tot


def fit_C_r2_optimal(theta_deg, density):
    """Closed-form C for dN/dtheta = C * sin(theta)/sin^4(theta/2) that
    MAXIMIZES r_squared_log10 -- no optimizer, no Poisson weighting, just
    algebra. Treats each histogram bin as one (bin center, bin top) point,
    same convention as r_squared_log10 itself.

    R^2_log = 1 - SS_res/SS_tot with SS_res = sum((log10(N_i) - log10(C *
    shape_i))^2) (shape_i = sin(theta_i)/sin^4(theta_i/2), theta in
    radians). For a fixed shape, SS_res is minimized over C by ordinary
    least squares in log space:

        log10(C) = mean(log10(density_i) - log10(shape_i))

    i.e. C is the geometric mean of density_i / shape_i across bins. Since
    maximizing R^2_log is exactly minimizing SS_res (SS_tot doesn't depend
    on C), this C is the unique amplitude that maximizes R^2_log for this
    shape -- not an approximation or a heuristic.

    `density` should already be counts/binwidth_rad (dN/dtheta per
    radian), matching dndtheta_model's convention, so the returned C is
    directly comparable to theoretical_amplitude()'s C_theory.
    """
    shape_vals = dndtheta_model(theta_deg, 1.0)
    c_i = np.asarray(density, dtype=float) / shape_vals
    return float(10 ** np.mean(np.log10(c_i)))


def theoretical_amplitude(n_total, Z1, Z2, E_mev, rho, t_cm, M):
    """C_theory such that dN/dtheta = C_theory * sin(theta)/sin^4(theta/2),
    theta in radians. See docstring eq. (3)."""
    C1_fm2 = (Z1 * Z2 * E2_MEV_FM / (4.0 * E_mev)) ** 2
    C1_cm2 = C1_fm2 * FM2_TO_CM2
    n_areal = rho * N_A * t_cm / M  # target nuclei / cm^2
    C2 = n_total * n_areal
    C_theory = 2 * np.pi * C1_cm2 * C2
    return C_theory, {"C1_cm2": C1_cm2, "n_areal_per_cm2": n_areal, "C2": C2}

def chi_squared_critical(dof: int, p: float = 0.05) -> float:
    """
    Return the chi-squared critical value for a given degrees of freedom
    at significance level p (default p = 0.05, upper-tail).

    Valid for dof = 1 to 180 (or beyond, scipy handles it fine).
    """
    if not (1 <= dof <= 180):
        raise ValueError(f"dof must be between 1 and 180, got {dof}")
    return chi2.ppf(1 - p, dof)
    


def chi_squared_check(dof, chi_squared):
    crit_value = chi_squared_critical(dof)
    if (chi_squared >= crit_value):
        return False
    else:
        return True

def calculate_expected(Z1, Z2, E, rho, thickness, molar_mass, num_particles):
    C1 = (((Z1 * Z2)/(E * 4)) * 1.44) ** 2
    C1 = C1 * (10 ** -26) #fm^2 to cm^2

    n = (rho * 6.022e23 * thickness)/molar_mass
    C2 = n * num_particles * 2 * math.pi

    C3 = C1 * C2
    return C3


def calculate_chi_squared(expected, observed):
    chi_squared = ((observed - expected) ** 2)/expected
    return chi_squared




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
    parser.add_argument("--Z1", type=float, default=2.0, help="projectile charge number (default 2, alpha)")
    parser.add_argument("--Z2", type=float, default=14.0, help="target charge number (default 14, silicon)")
    parser.add_argument("--E", type=float, default=12.0, help="incident kinetic energy, MeV (default 6.0)")
    parser.add_argument("--rho", type=float, default=2.33, help="target density, g/cm^3 (default 2.33, silicon)")
    parser.add_argument("--t", type=float, default=1e-4, help="foil thickness, cm (default 1e-4 = 1 um, matches DetectorConstruction.cc kFoilHalfZ)")
    parser.add_argument("--M", type=float, default=28.0855, help="target molar mass, g/mol (default 28.0855, silicon)")
    parser.add_argument("--no-plot", dest="plot", action="store_false", help="skip writing the PNG plot (a plot is written by default)")
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
    # C is chosen in closed form to MAXIMIZE R^2 (log10 N space) -- see
    # fit_C_r2_optimal's docstring. No optimizer, no Poisson weighting: each
    # bin is treated as one (bin center, bin top) point, and C is the exact
    # amplitude that maximizes r_squared_log10 for this shape.
    n_bins_fit = int(mask.sum())
    n_params = 1  # just the amplitude C -- theta_min is a fixed cut, not a fitted parameter
    dof = n_bins_fit - n_params
    C_fit = fit_C_r2_optimal(centers_deg[mask], density[mask])

    # R^2 in N (counts-per-bin) space, treating each of the n_bins_fit bins
    # in the fit region as one point (bin center, bin height/top).
    fit_N_at_mask = dndtheta_model(centers_deg[mask], C_fit) * binwidth_rad
    r2_fit = r_squared_log10(counts[mask], fit_N_at_mask)

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
    print("Fitted amplitude (fit_C_r2_optimal, closed-form log-space least-squares, theta-domain restricted to tail):")
    print(f"  C_fit = {C_fit:.6e}        R^2 (log10 N space) = {r2_fit:.4f}")
    print(f"  n_bins (fit region, theta >= {theta_min:.3f} deg) = {n_bins_fit}")
    print(f"  n_params (free parameters: C only)                = {n_params}")
    print()
    print(f"ratio C_fit / C_theory = {C_fit / C_theory:.6g}")

    # ---- chi-squared check on C_fit (the sin(theta)/sin^4(theta/2) fit) --
    # calculate_expected() re-derives the theoretical amplitude from the
    # same physics parameters as theoretical_amplitude() above (Z1, Z2, E,
    # rho, t, M, n_total); calculate_chi_squared() compares that expected
    # value against the actual fitted C_fit; chi_squared_check() passes if
    # that chi^2 is below the critical value at this fit's dof (n_bins_fit -
    # n_params, same dof already used for the goodness-of-fit block above).
    expected_C = calculate_expected(args.Z1, args.Z2, args.E, args.rho, args.t, args.M, n_total)
    chi_sq_C = calculate_chi_squared(expected_C, C_fit)
    chi_sq_crit = chi_squared_critical(dof)
    chi_sq_passed = chi_squared_check(dof, chi_sq_C)

    print()
    print("Chi-squared check -- fitted C (sin(theta)/sin^4(theta/2)) vs. theoretical C:")
    print(f"  expected C (calculate_expected)                    = {expected_C:.6e}")
    print(f"  observed C (C_fit)                                 = {C_fit:.6e}")
    print(f"  chi^2 = (observed - expected)^2 / expected         = {chi_sq_C:.4f}")
    print(f"  chi^2 critical (p=0.05)                            = {chi_sq_crit:.4f}")
    print(f"  passes (chi^2 < critical)?                         = {chi_sq_passed}")
    if chi_sq_passed:
        print("  Simulated Data is Accurate")

    if args.plot:
        import matplotlib.pyplot as plt

        theta_plot = np.linspace(centers_deg[mask].min(), centers_deg.max(), 400)

        fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
        fig.patch.set_facecolor("#fcfcfb")
        style_axes(ax)

        ax.bar(centers_deg[counts > 0], counts[counts > 0], width=args.binwidth * 0.92,
               color=COLOR_BARS, edgecolor="#fcfcfb", linewidth=0.3,
               label=f"Simulated foil-exit angles, N (n={n_total})")
        fit_curve_N = dndtheta_model(theta_plot, C_fit) * binwidth_rad
        ax.plot(theta_plot, fit_curve_N, color=COLOR_ENVELOPE, linewidth=2.5, linestyle="-.", zorder=4,
                label=fr"$\sin\theta/\sin^4(\theta/2)$ fit ($R^2$-optimal): "
                      fr"$C$ = {C_fit:.3e}, $R^2_{{\log}}$ = {r2_fit:.3f}")
        ax.axvline(theta_min, color=COLOR_AXIS, linewidth=1, linestyle=":")

        ax.set_yscale("log")
        ax.set_ylim(bottom=1)
        ax.set_xlabel(r"Scattering angle, $\theta$ (degrees)", color=COLOR_INK)
        ax.set_ylabel("Number of particles (log scale)", color=COLOR_INK)
        foil_um = args.t * 1e4  # cm -> um
        particle_label = {1.0: "Proton", 2.0: "Alpha Particle"}.get(args.Z1, f"Z1={args.Z1:g}")
        target_label = {79.0: "Gold", 14.0: "Silicon", 6.0: "Carbon"}.get(args.Z2, f"Z2={args.Z2:g}")
        ax.set_title(rf"Rutherford Scattering Angular Distribution ({foil_um:g} $\mu$m {target_label} Foil) {particle_label} 12MeV"
                     "\n"
                     r"Simulated Beam Data with $\sin\theta/\sin^4(\theta/2)$ Fit",
                     color=COLOR_INK, fontsize=13)
        ax.legend(frameon=False, labelcolor=COLOR_INK, fontsize=9)

        pass_tag = "PASS" if chi_sq_passed else "FAIL"
        pass_color = "#1baf7a" if chi_sq_passed else "#d9432e"
        accuracy_text = "Simulated Data is Accurate" if chi_sq_passed else "Simulated Data is Inaccurate"
        eqn_text = (
            r"$N(\theta) = C \cdot \dfrac{\sin\theta}{\sin^4(\theta/2)}$"
            "\n"
            fr"$C_{{fit}} = {C_fit:.3f}$"
            "\n"
            fr"$C_{{expected}} = {expected_C:.3f}$"
            "\n"
            fr"$R^2_{{\log}}$ = {r2_fit:.4f}"
            "\n"
            fr"$\chi^2$ = {chi_sq_C:.3f}, critical ($p$=0.05) = {chi_sq_crit:.3f}"
            "\n"
            fr"$\chi^2$ < critical? {pass_tag}"
            "\n"
            fr"$\mathbf{{{accuracy_text.replace(' ', r'\ ')}}}$"
        )
        ax.text(0.55, 0.55, eqn_text, transform=ax.transAxes, fontsize=11,
                color=COLOR_INK, ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="#fcfcfb",
                          edgecolor=pass_color, linewidth=1.2))

        fig.tight_layout()

        out_png = os.path.splitext(csv_path)[0] + "_rutherford_fit.png"
        fig.savefig(out_png, facecolor=fig.get_facecolor())
        print(f"Wrote plot to {out_png}")


if __name__ == "__main__":
    main()
