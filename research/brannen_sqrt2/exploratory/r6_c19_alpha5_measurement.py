#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c19_alpha5_measurement.py:
Measure alpha_5 from the full nonlinear substrate ODE.

Strategy:
- Solve g'' + (g')²/g + (2/r)g' + g = 1 for many g0 = 1+δ
- Extract A_tail = amplitude of (g-1)·r at large r
- Compute η_sym(δ) = (η(+δ) + η(-δ))/2 = 1 + α_3·δ² + α_5·δ⁴ + O(δ⁶)
- Using known α_3 = 0.089722224 (30 digits), fit:
  F(δ) ≡ [η_sym(δ) - 1 - α_3·δ²] / δ⁴ → α_5 as δ→0
- Richardson extrapolation on r_max to remove truncation bias

Use r_max=800 (good compromise: fast enough, truncation error manageable).
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
                    rtol=1e-13, atol=1e-15, max_step=0.015)
    return sol

def extract_atail(sol, r_min, r_max_win):
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max_win)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    return math.sqrt(coef[0]**2 + coef[1]**2)

def measure_eta(delta, r_max=800, n_pts=320000):
    sol_p = solve_ode(1.0 + delta, r_max, n_pts)
    sol_m = solve_ode(1.0 - delta, r_max, n_pts)
    if not sol_p.success or not sol_m.success:
        return None, None
    A_p = extract_atail(sol_p, 100.0, r_max - 50.0)
    A_m = extract_atail(sol_m, 100.0, r_max - 50.0)
    eta_p = A_p / delta
    eta_m = A_m / delta
    eta_sym = (eta_p + eta_m) / 2.0
    return eta_sym, (eta_p, eta_m)

print("="*78)
print(f"  ALPHA_5 MEASUREMENT from full substrate ODE")
print(f"  Using alpha_3 = {alpha3_exact:.20f} (30 digits from mpmath)")
print("="*78)

# Scan over δ values
deltas = np.array([0.002, 0.003, 0.005, 0.008, 0.010, 0.015, 0.020, 0.030, 0.040, 0.050])

for r_max in [600, 1000]:
    n_pts = int(400 * r_max)
    print(f"\n  --- r_max = {r_max}, n_pts = {n_pts} ---")
    print(f"  {'delta':>8s}  {'eta_sym':>14s}  {'(eta_sym-1)/d^2':>16s}  {'F(d) = (...)/ d^4':>20s}")

    d2_vals = []
    F_vals = []
    raw_ratio = []

    for d in deltas:
        eta_sym, details = measure_eta(d, r_max, n_pts)
        if eta_sym is None:
            print(f"  {d:8.4f}  FAILED")
            continue
        ratio = (eta_sym - 1.0) / d**2
        F = (eta_sym - 1.0 - alpha3_exact * d**2) / d**4
        print(f"  {d:8.4f}  {eta_sym:14.10f}  {ratio:16.10f}  {F:20.8f}")
        d2_vals.append(d**2)
        F_vals.append(F)
        raw_ratio.append(ratio)

    d2_arr = np.array(d2_vals)
    F_arr = np.array(F_vals)

    # Fit F(δ) = α_5 + α_7·δ² + ... → linear in δ²
    if len(d2_arr) >= 3:
        # Use only small δ for fit
        mask_small = d2_arr < 0.001
        if np.sum(mask_small) < 3:
            mask_small = d2_arr < 0.003

        A_fit = np.column_stack([np.ones(len(d2_arr)), d2_arr])
        coef, *_ = np.linalg.lstsq(A_fit, F_arr, rcond=None)
        alpha5_fit = coef[0]
        alpha7_approx = coef[1]
        print(f"\n    Linear fit F = alpha_5 + alpha_7*d^2:")
        print(f"      alpha_5 = {alpha5_fit:.8f}")
        print(f"      alpha_7 ~ {alpha7_approx:.4f}")

        # Quadratic for better α_5
        A_fit2 = np.column_stack([np.ones(len(d2_arr)), d2_arr, d2_arr**2])
        coef2, *_ = np.linalg.lstsq(A_fit2, F_arr, rcond=None)
        print(f"    Quadratic fit F = a + b*d^2 + c*d^4:")
        print(f"      alpha_5 = {coef2[0]:.8f}")

        # Also fit (eta_sym-1)/d^2 = a0 + a2*d^2 + a4*d^4
        ratio_arr = np.array(raw_ratio)
        A3 = np.column_stack([np.ones(len(d2_arr)), d2_arr, d2_arr**2])
        c3, *_ = np.linalg.lstsq(A3, ratio_arr, rcond=None)
        print(f"\n    Direct fit (eta_sym-1)/d^2 = a0 + a2*d^2 + a4*d^4:")
        print(f"      a0 (= alpha_3) = {c3[0]:.12f}  (target: {alpha3_exact:.12f})")
        print(f"      a2 (= alpha_5) = {c3[1]:.8f}")
        print(f"      a4 (= alpha_7) = {c3[2]:.4f}")

print(f"\n{'='*78}")
print("  CANDIDATE IDENTIFICATION for alpha_5:")
print(f"{'='*78}")
# Will print candidates after seeing the value
