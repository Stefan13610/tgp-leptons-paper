#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c26_alpha5_large_rmax.py:
Alpha_5 from full ODE with very large R_max (2000, 3000) for minimal truncation.

Key improvements:
1. Much larger R_max → truncation error shrinks as 1/R_max
2. Large fitting window → better amplitude extraction
3. Include 1/r decay correction in fit
4. Use many delta values (0.06-0.25) for robust polynomial fit
5. Multiple R_max for Richardson extrapolation
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
import math
from scipy.integrate import solve_ivp

alpha3_exact = 0.08972222367362532604749

def rhs(r, y):
    g, gp = y
    if g < 1e-12: g = 1e-12
    if r < 1e-13:
        gpp = (1.0 - g) / 4.0
    else:
        gpp = (1.0 - g) - (1.0/g) * gp**2 - (2.0/r) * gp
    return [gp, gpp]

def solve_and_extract(g0, r_max, n_pts, fit_lo, fit_hi):
    """Solve ODE and extract tail amplitude."""
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=2.3e-14, atol=1e-15, max_step=0.05)
    if not sol.success:
        return None
    r = sol.t
    g = sol.y[0]
    mask = (r >= fit_lo) & (r <= fit_hi)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    # Fit: u = Ac*cos(r) + As*sin(r) + Bc*cos(r)/r + Bs*sin(r)/r
    X = np.column_stack([
        np.cos(r_f), np.sin(r_f),
        np.cos(r_f)/r_f, np.sin(r_f)/r_f
    ])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    A = math.sqrt(coef[0]**2 + coef[1]**2)
    return A

print("="*78)
print("  ALPHA_5 with very large R_max (2000, 3000)")
print("="*78)

deltas = np.array([0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25])

configs = [
    # (r_max, n_pts, fit_lo, fit_hi)
    (2000, 600000, 500, 1900),
    (3000, 800000, 800, 2800),
]

all_results = {}

for r_max, n_pts, fit_lo, fit_hi in configs:
    print(f"\n  --- R_max={r_max}, fit=[{fit_lo},{fit_hi}], n_pts={n_pts} ---")
    print(f"  {'d':>6s} {'eta_sym':>18s} {'ratio':>18s} {'F':>16s}")

    d2_list = []
    ratio_list = []

    for d in deltas:
        Ap = solve_and_extract(1.0+d, r_max, n_pts, fit_lo, fit_hi)
        Am = solve_and_extract(1.0-d, r_max, n_pts, fit_lo, fit_hi)
        if Ap is None or Am is None:
            print(f"  {d:6.2f} FAILED")
            continue
        eta_p = Ap / d
        eta_m = Am / d
        eta_sym = (eta_p + eta_m) / 2.0
        ratio = (eta_sym - 1.0) / d**2
        F = (eta_sym - 1.0 - alpha3_exact * d**2) / d**4
        print(f"  {d:6.2f} {eta_sym:18.14f} {ratio:18.14f} {F:16.10f}")
        d2_list.append(d**2)
        ratio_list.append(ratio)

    d2 = np.array(d2_list)
    ratio = np.array(ratio_list)

    if len(d2) < 5:
        continue

    # Polynomial fits of ratio = alpha_3 + alpha_5*d^2 + alpha_7*d^4 + ...
    for deg in [2, 3, 4]:
        if len(d2) < deg + 2:
            continue
        A_mat = np.column_stack([d2**k for k in range(deg+1)])
        coef, res, *_ = np.linalg.lstsq(A_mat, ratio, rcond=None)

        labels = ['alpha_3', 'alpha_5', 'alpha_7', 'alpha_9', 'alpha_11']
        vals = []
        for i in range(min(deg+1, 5)):
            vals.append(f"{labels[i]}={coef[i]:.14f}")
        a3_err = coef[0] - alpha3_exact
        print(f"    deg {deg}: {', '.join(vals)}  (a3 err={a3_err:.3e})")

    all_results[r_max] = {
        'd2': d2.copy(),
        'ratio': ratio.copy(),
    }

# Richardson extrapolation on alpha_5
print(f"\n{'='*78}")
print(f"  RICHARDSON EXTRAPOLATION")
print(f"{'='*78}")

rmax_list = sorted(all_results.keys())
if len(rmax_list) >= 2:
    # Use degree 3 polynomial fit for each R_max
    a5_by_rmax = {}
    for rm in rmax_list:
        d2 = all_results[rm]['d2']
        ratio = all_results[rm]['ratio']
        A_mat = np.column_stack([d2**k for k in range(4)])
        coef, *_ = np.linalg.lstsq(A_mat, ratio, rcond=None)
        a5_by_rmax[rm] = coef[1]
        print(f"  R_max={rm}: alpha_5 = {coef[1]:.14f}, alpha_3_err = {coef[0]-alpha3_exact:.3e}")

    # Richardson (1/R model)
    r1, r2 = rmax_list
    a1, a2 = a5_by_rmax[r1], a5_by_rmax[r2]
    a5_rich_1r = (r2*a2 - r1*a1) / (r2 - r1)
    print(f"\n  Richardson 1/R: alpha_5 = {a5_rich_1r:.14f}")

    # Richardson (1/R^2 model)
    a5_rich_2r = (r2**2*a2 - r1**2*a1) / (r2**2 - r1**2)
    print(f"  Richardson 1/R^2: alpha_5 = {a5_rich_2r:.14f}")

# Summary
print(f"\n{'='*78}")
print(f"  SUMMARY")
print(f"{'='*78}")
print(f"  Previous estimates: scipy R_max=400-1000 → alpha_5 ≈ 0.0275 (4 digits)")
print(f"  This run: large R_max → reduced truncation → better convergence")
print(f"  Best estimate from Richardson: see above")
print(f"  Known: alpha_3 = {alpha3_exact:.20f}")
print("="*78)
