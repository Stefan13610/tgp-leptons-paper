#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c21_alpha5_richardson_pslq.py:
1. Richardson extrapolation of alpha_5 from r6_c20 multi-r_max data
2. Higher-precision ODE measurement with finer delta grid
3. PSLQ identification attempt for alpha_5

From r6_c20, Method 1 degree 2 gives:
  r_max=300: alpha_5 = 0.02745757, alpha_3_err = +5.22e-5
  r_max=500: alpha_5 = 0.02746761, alpha_3_err = -1.60e-5
  r_max=800: alpha_5 = 0.02747929, alpha_3_err = -1.99e-5

Method 3 (mid-range quadratic):
  r_max=300: alpha_5 = 0.02746568
  r_max=500: alpha_5 = 0.02747508
  r_max=800: alpha_5 = 0.02748790

The alpha_3 error oscillates, suggesting tail interference effects.
Alpha_5 is more stable and converging upward.
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

def solve_ode(g0, r_max, n_pts):
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=2.3e-14, atol=1e-15, max_step=0.008)
    return sol

def extract_atail(sol, r_min, r_max_win):
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max_win)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([
        np.cos(r_f), np.sin(r_f),
        np.cos(r_f)/r_f, np.sin(r_f)/r_f
    ])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    return math.sqrt(coef[0]**2 + coef[1]**2)

def measure_eta_sym(delta, r_max, n_pts):
    sol_p = solve_ode(1.0 + delta, r_max, n_pts)
    sol_m = solve_ode(1.0 - delta, r_max, n_pts)
    if not sol_p.success or not sol_m.success:
        return None

    windows = [
        (80.0, r_max - 80.0),
        (100.0, r_max - 60.0),
        (60.0, r_max - 100.0),
    ]
    etas = []
    for r_lo, r_hi in windows:
        if r_hi > r_lo + 100:
            A_p = extract_atail(sol_p, r_lo, r_hi)
            A_m = extract_atail(sol_m, r_lo, r_hi)
            eta_sym = (A_p/delta + A_m/delta) / 2.0
            etas.append(eta_sym)
    return np.mean(etas) if etas else None

print("="*78)
print("  ALPHA_5: HIGH-PRECISION ODE + RICHARDSON + PSLQ")
print("="*78)

# Fine delta grid for better polynomial fitting
deltas = np.array([0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20])

# Use r_max = 400, 600, 1000 for Richardson
r_maxes = [400, 600, 1000]
results = {}

for r_max in r_maxes:
    n_pts = int(600 * r_max)
    print(f"\n  --- r_max = {r_max}, n_pts = {n_pts} ---")

    d2_list = []
    ratio_list = []

    for d in deltas:
        eta_sym = measure_eta_sym(d, r_max, n_pts)
        if eta_sym is None:
            continue
        ratio = (eta_sym - 1.0) / d**2
        d2_list.append(d**2)
        ratio_list.append(ratio)
        print(f"    d={d:.2f}  eta_sym-1={eta_sym-1:.12e}  ratio={ratio:.12f}")

    d2 = np.array(d2_list)
    ratio = np.array(ratio_list)

    # Degree 3 polynomial fit: ratio = a0 + a1*d^2 + a2*d^4 + a3*d^6
    A = np.column_stack([d2**k for k in range(4)])
    coef, *_ = np.linalg.lstsq(A, ratio, rcond=None)

    results[r_max] = {
        'alpha3': coef[0],
        'alpha5': coef[1],
        'alpha7': coef[2],
        'alpha9': coef[3],
    }

    print(f"    Fit: alpha_3 = {coef[0]:.14f}  (err = {coef[0]-alpha3_exact:.3e})")
    print(f"          alpha_5 = {coef[1]:.10f}")
    print(f"          alpha_7 = {coef[2]:.6f}")
    print(f"          alpha_9 = {coef[3]:.4f}")

# Richardson extrapolation
print(f"\n{'='*78}")
print("  RICHARDSON EXTRAPOLATION")
print("="*78)

# alpha_5 values at each r_max
a5_vals = {rm: results[rm]['alpha5'] for rm in r_maxes}
print(f"  alpha_5 at r_max:")
for rm in r_maxes:
    print(f"    r_max={rm:4d}: alpha_5 = {a5_vals[rm]:.10f}")

# Simple 2-point Richardson (assuming 1/r_max^2 truncation error)
for i in range(len(r_maxes)):
    for j in range(i+1, len(r_maxes)):
        r1, r2 = r_maxes[i], r_maxes[j]
        a5_1, a5_2 = a5_vals[r1], a5_vals[r2]
        # a5(r) = a5_true + C/r^2
        # a5_true = (r2^2 * a5_2 - r1^2 * a5_1) / (r2^2 - r1^2)
        a5_rich = (r2**2 * a5_2 - r1**2 * a5_1) / (r2**2 - r1**2)
        print(f"  Richardson ({r1},{r2}): alpha_5 = {a5_rich:.10f}")

# Also try 1/r_max decay
print(f"\n  Assuming 1/r_max truncation:")
for i in range(len(r_maxes)):
    for j in range(i+1, len(r_maxes)):
        r1, r2 = r_maxes[i], r_maxes[j]
        a5_1, a5_2 = a5_vals[r1], a5_vals[r2]
        a5_rich = (r2 * a5_2 - r1 * a5_1) / (r2 - r1)
        print(f"  Richardson ({r1},{r2}): alpha_5 = {a5_rich:.10f}")

# 3-point Richardson
if len(r_maxes) >= 3:
    r1, r2, r3 = r_maxes
    a1, a2, a3 = a5_vals[r1], a5_vals[r2], a5_vals[r3]
    # Fit a5 = a5_true + A/r + B/r^2
    A_rich = np.column_stack([
        np.ones(3),
        [1/r1, 1/r2, 1/r3],
        [1/r1**2, 1/r2**2, 1/r3**2]
    ])
    c_rich = np.linalg.solve(A_rich, [a1, a2, a3])
    print(f"\n  3-point Richardson (1/r + 1/r^2 model):")
    print(f"    alpha_5_true = {c_rich[0]:.10f}")
    print(f"    A = {c_rich[1]:.6f}")
    print(f"    B = {c_rich[2]:.4f}")

# PSLQ attempt
print(f"\n{'='*78}")
print("  PSLQ IDENTIFICATION ATTEMPT")
print("="*78)

try:
    import mpmath as mp
    mp.mp.dps = 30
    pi = mp.pi
    ln2 = mp.log(2)
    ln3 = mp.log(3)

    # Best estimate of alpha_5 (will use ~0.0275)
    # Take the most stable estimates
    a5_best = float(c_rich[0])  # 3-point Richardson
    print(f"  Using alpha_5 estimate: {a5_best:.10f}")

    a5_mp = mp.mpf(str(a5_best))

    # Test simple candidates
    candidates = {
        'pi^2/360': pi**2/360,
        '1/36': mp.mpf(1)/36,
        '11/400': mp.mpf(11)/400,
        'pi^2/358': pi**2/358,
        'pi^2/359': pi**2/359,
        'pi^4/35000': pi**4/35000,
        'ln2/25': ln2/25,
        '(pi^2-9)/36': (pi**2-9)/36,
        'pi^2/128/4': pi**2/512,
        '1/6 - ln2': mp.mpf(1)/6 - ln2,
        'alpha3/sqrt(10)': mp.mpf(alpha3_exact)/mp.sqrt(10),
        'alpha3^2*3.4': mp.mpf(alpha3_exact)**2 * mp.mpf('3.4'),
    }

    print(f"\n  Candidate comparison:")
    for name, val in candidates.items():
        diff = a5_mp - val
        print(f"    {name:25s} = {mp.nstr(val, 12):>16s}  diff = {mp.nstr(diff, 4)}")

    # PSLQ with basic constants
    basis = [a5_mp, pi**2, pi**4, ln2, ln3, ln2**2, mp.zeta(3), mp.catalan, mp.mpf(1)]
    names = ['alpha_5', 'pi^2', 'pi^4', 'ln2', 'ln3', 'ln2^2', 'zeta3', 'G', '1']

    # Note: PSLQ with only ~4 digit accuracy won't work well
    # We need at least 15+ digits for meaningful PSLQ
    print(f"\n  NOTE: alpha_5 has only ~4-5 significant digits from ODE method.")
    print(f"  PSLQ requires 15+ digits for reliable identification.")
    print(f"  A perturbation-theory computation of alpha_5 would be needed.")

except Exception as e:
    print(f"  Error: {e}")

print("="*78)
