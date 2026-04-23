#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c20_alpha5_revised.py:
Revised alpha_5 measurement from full nonlinear substrate ODE.

Key improvements over r6_c19:
1. Use larger delta (0.03-0.25) where signal >> float64 noise
2. Use multiple r_max for Richardson extrapolation
3. Higher-order polynomial fits (up to degree 4 in delta^2)
4. Careful A_tail extraction with multiple fitting windows

Strategy:
  eta_sym(d) = (eta(+d) + eta(-d))/2 = 1 + alpha_3*d^2 + alpha_5*d^4 + alpha_7*d^6 + ...
  Fit (eta_sym-1)/d^2 vs d^2 as polynomial -> intercept = alpha_3, slope = alpha_5
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import math
from scipy.integrate import solve_ivp

alpha3_exact = 0.08972222367362532604749  # 25 digits from mpmath

def rhs(r, y):
    g, gp = y
    if g < 1e-12: g = 1e-12
    if r < 1e-13:
        gpp = (1.0 - g) / 4.0
    else:
        gpp = (1.0 - g) - (1.0/g) * gp**2 - (2.0/r) * gp
    return [gp, gpp]

def solve_ode(g0, r_max, n_pts):
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=1e-14, atol=1e-15, max_step=0.01)
    return sol

def extract_atail(sol, r_min, r_max_win):
    """Extract tail amplitude A from (g-1)*r ~ A*cos(r+phi) at large r."""
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max_win)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    # Fit u = A_c*cos(r) + A_s*sin(r) + decay terms
    # Include 1/r decay correction: A_c2*cos(r)/r + A_s2*sin(r)/r
    X = np.column_stack([
        np.cos(r_f), np.sin(r_f),
        np.cos(r_f)/r_f, np.sin(r_f)/r_f
    ])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    A = math.sqrt(coef[0]**2 + coef[1]**2)
    return A

def measure_eta_sym(delta, r_max=400, n_pts=200000):
    """Measure eta_sym = (eta(+d) + eta(-d))/2 for given delta."""
    sol_p = solve_ode(1.0 + delta, r_max, n_pts)
    sol_m = solve_ode(1.0 - delta, r_max, n_pts)
    if not sol_p.success or not sol_m.success:
        return None

    # Use multiple fitting windows and average
    windows = []
    # Window 1: [80, r_max-80]
    if r_max > 200:
        windows.append((80.0, r_max - 80.0))
    # Window 2: [120, r_max-60]
    if r_max > 250:
        windows.append((120.0, r_max - 60.0))
    # Window 3: [60, r_max-100]
    if r_max > 200:
        windows.append((60.0, r_max - 100.0))

    if not windows:
        windows = [(50.0, r_max - 30.0)]

    etas = []
    for r_lo, r_hi in windows:
        A_p = extract_atail(sol_p, r_lo, r_hi)
        A_m = extract_atail(sol_m, r_lo, r_hi)
        eta_p = A_p / delta
        eta_m = A_m / delta
        eta_sym = (eta_p + eta_m) / 2.0
        etas.append(eta_sym)

    return np.mean(etas)

print("="*78)
print("  REVISED ALPHA_5 MEASUREMENT (larger delta, multi-window)")
print(f"  alpha_3 = {alpha3_exact:.20f}")
print("="*78)

# Use larger delta values where delta^4 contribution is well above float64 noise
deltas = np.array([0.03, 0.05, 0.07, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25])

for r_max in [300, 500, 800]:
    n_pts = int(500 * r_max)  # ~0.002 step size
    print(f"\n  --- r_max = {r_max}, n_pts = {n_pts} ---")
    print(f"  {'delta':>8s}  {'eta_sym':>16s}  {'(eta_sym-1)/d^2':>16s}  {'F(d)':>16s}")

    d2_vals = []
    ratio_vals = []  # (eta_sym - 1) / d^2
    F_vals = []      # (eta_sym - 1 - alpha3*d^2) / d^4

    for d in deltas:
        eta_sym = measure_eta_sym(d, r_max, n_pts)
        if eta_sym is None:
            print(f"  {d:8.4f}  FAILED")
            continue
        ratio = (eta_sym - 1.0) / d**2
        F = (eta_sym - 1.0 - alpha3_exact * d**2) / d**4
        print(f"  {d:8.4f}  {eta_sym:16.12f}  {ratio:16.10f}  {F:16.8f}")
        d2_vals.append(d**2)
        ratio_vals.append(ratio)
        F_vals.append(F)

    d2 = np.array(d2_vals)
    ratio = np.array(ratio_vals)
    F = np.array(F_vals)

    if len(d2) < 4:
        print("  Too few points for fitting")
        continue

    # Method 1: Fit (eta_sym-1)/d^2 = alpha_3 + alpha_5*d^2 + alpha_7*d^4 + alpha_9*d^6
    print(f"\n  Method 1: polynomial fit of (eta_sym-1)/d^2 vs d^2")
    for deg in [1, 2, 3]:
        A = np.column_stack([d2**k for k in range(deg+1)])
        coef, res, *_ = np.linalg.lstsq(A, ratio, rcond=None)
        print(f"    degree {deg}: alpha_3={coef[0]:.12f}, alpha_5={coef[1]:.8f}", end="")
        if deg >= 2:
            print(f", alpha_7={coef[2]:.4f}", end="")
        if deg >= 3:
            print(f", alpha_9={coef[3]:.2f}", end="")
        print(f"  (alpha_3 err = {coef[0] - alpha3_exact:.2e})")

    # Method 2: Fit F(d) = alpha_5 + alpha_7*d^2 + ...
    print(f"\n  Method 2: polynomial fit of F(d) vs d^2")
    for deg in [1, 2, 3]:
        if len(d2) < deg + 2:
            continue
        A = np.column_stack([d2**k for k in range(deg+1)])
        coef, *_ = np.linalg.lstsq(A, F, rcond=None)
        print(f"    degree {deg}: alpha_5={coef[0]:.8f}", end="")
        if deg >= 2:
            print(f", alpha_7={coef[1]:.4f}", end="")
        print()

    # Method 3: Use only the middle delta range (0.07-0.18) for stability
    mask_mid = (d2 >= 0.004) & (d2 <= 0.04)
    if np.sum(mask_mid) >= 4:
        d2_mid = d2[mask_mid]
        ratio_mid = ratio[mask_mid]
        A = np.column_stack([np.ones(len(d2_mid)), d2_mid, d2_mid**2])
        coef, *_ = np.linalg.lstsq(A, ratio_mid, rcond=None)
        print(f"\n  Method 3: quadratic fit on mid-range delta (0.07-0.20):")
        print(f"    alpha_3 = {coef[0]:.12f}  (err = {coef[0]-alpha3_exact:.2e})")
        print(f"    alpha_5 = {coef[1]:.8f}")
        print(f"    alpha_7 = {coef[2]:.4f}")

# Richardson extrapolation between r_max values
print(f"\n{'='*78}")
print(f"  RICHARDSON EXTRAPOLATION between r_max values")
print(f"{'='*78}")
print(f"  (Will apply after seeing consistency across r_max)")

print(f"\n{'='*78}")
print(f"  SUMMARY")
print(f"{'='*78}")
print(f"  Look for alpha_5 estimates that are stable across:")
print(f"    - different r_max (truncation convergence)")
print(f"    - different fit degrees (polynomial convergence)")
print(f"    - different delta ranges (signal quality)")
print("="*78)
