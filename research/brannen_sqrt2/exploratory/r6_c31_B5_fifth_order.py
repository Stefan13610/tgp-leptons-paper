#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c31_B5_fifth_order.py:
Compute B5, A5 (and reconfirm B4, A4) from 5th-order perturbation ODE.

From g'' + (1/g)(g')^2 + (2/r)g' + g = 1 with g = 1 + h, h = sum d^n f_n:
  h'' + (h')^2/(1+h) + (2/r)h' + h = 0
  Expanding 1/(1+h) = 1 - h + h^2 - h^3 + h^4 - ...

Source terms (f_n'' + (2/r)f_n' + f_n = S_n):

S2 = -(f1')^2

S3 = f1*(f1')^2 - 2*f1'*f2'

S4 = 2*f1*f1'*f2' + f2*(f1')^2 - f1^2*(f1')^2 - 2*f1'*f3' - (f2')^2

S5 = -2*f1'*f4' - 2*f2'*f3'
     + 2*f1*f1'*f3' + f1*(f2')^2 + 2*f2*f1'*f2' + f3*(f1')^2
     - 2*f1^2*f1'*f2' - 2*f1*f2*(f1')^2
     + f1^3*(f1')^2

(Derivation: collect delta^5 terms from (h')^2 * (1 - h + h^2 - h^3 + ...))

ODE system: 8 components [f2, f2', f3, f3', f4, f4', f5, f5']
f1 = sin(r)/r is known analytically.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

print("="*78)
print("  B5, A5 via 5th-order perturbation ODE (scipy DOP853)")
print("="*78)

def f1(r):
    if r < 1e-10:
        return 1.0 - r**2/6 + r**4/120 - r**6/5040
    return np.sin(r)/r

def f1p(r):
    if r < 1e-10:
        return -r/3 + r**3/30 - r**5/840 + r**7/45360
    return np.cos(r)/r - np.sin(r)/r**2

def rhs(r, y):
    f2, f2p, f3, f3p, f4, f4p, f5, f5p = y

    if r < 1e-10:
        return [f2p, 0, f3p, 0, f4p, 0, f5p, 0]

    fv = f1(r)
    fp = f1p(r)
    tr = 2.0/r

    # Source terms
    s2 = -fp**2
    s3 = fv*fp**2 - 2*fp*f2p
    s4 = (2*fv*fp*f2p + f2*fp**2 - fv**2*fp**2
          - 2*fp*f3p - f2p**2)
    s5 = (-2*fp*f4p - 2*f2p*f3p
          + 2*fv*fp*f3p + fv*f2p**2 + 2*f2*fp*f2p + f3*fp**2
          - 2*fv**2*fp*f2p - 2*fv*f2*fp**2
          + fv**3*fp**2)

    return [f2p, s2 - tr*f2p - f2,
            f3p, s3 - tr*f3p - f3,
            f4p, s4 - tr*f4p - f4,
            f5p, s5 - tr*f5p - f5]

r0 = 1e-6
y0 = [0.0]*8

import math
ln3 = math.log(3)
B2_known = 0.5 - ln3/8
A2_known = math.pi/8
B3_known = 0.012615939290114712
A3_known = -0.215712018451450166

def extract_amplitudes(r_arr, u_arr, r_fit_min, r_fit_max):
    mask = (r_arr >= r_fit_min) & (r_arr <= r_fit_max)
    r = r_arr[mask]
    u = u_arr[mask]
    M = np.column_stack([np.sin(r), np.cos(r), np.sin(r)/r, np.cos(r)/r])
    coeffs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    return coeffs[0], coeffs[1]

# Multiple R_max for Richardson
R_max_list = [1000, 2000, 3000]
results = {}

for R_max in R_max_list:
    print(f"\n  --- R_max = {R_max} ---")
    r_eval = np.linspace(r0, R_max, int(R_max * 20) + 1)

    sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval,
                    max_step=0.1)

    if not sol.success:
        print(f"    Failed: {sol.message}")
        continue

    r_arr = sol.t
    u2 = r_arr * sol.y[0]
    u3 = r_arr * sol.y[2]
    u4 = r_arr * sol.y[4]
    u5 = r_arr * sol.y[6]

    r_fit_min = R_max * 0.4
    r_fit_max = R_max * 0.95

    B2, A2 = extract_amplitudes(r_arr, u2, r_fit_min, r_fit_max)
    B3, A3 = extract_amplitudes(r_arr, u3, r_fit_min, r_fit_max)
    B4, A4 = extract_amplitudes(r_arr, u4, r_fit_min, r_fit_max)
    B5, A5 = extract_amplitudes(r_arr, u5, r_fit_min, r_fit_max)

    print(f"    B2 = {B2:+.14f}  (err: {B2-B2_known:+.2e})")
    print(f"    A2 = {A2:+.14f}  (err: {A2-A2_known:+.2e})")
    print(f"    B3 = {B3:+.14f}  (err: {B3-B3_known:+.2e})")
    print(f"    A3 = {A3:+.14f}  (err: {A3-A3_known:+.2e})")
    print(f"    B4 = {B4:+.14f}")
    print(f"    A4 = {A4:+.14f}")
    print(f"    B5 = {B5:+.14f}")
    print(f"    A5 = {A5:+.14f}")

    results[R_max] = {'B2':B2, 'A2':A2, 'B3':B3, 'A3':A3,
                      'B4':B4, 'A4':A4, 'B5':B5, 'A5':A5}

# Richardson
print(f"\n{'='*78}")
print("  RICHARDSON EXTRAPOLATION")
print("="*78)

Rvals = sorted(results.keys())
for amp in ['B4', 'A4', 'B5', 'A5']:
    print(f"\n  {amp}:")
    for R in Rvals:
        print(f"    R={R:5d}: {results[R][amp]:+.14f}")
    if len(Rvals) >= 2:
        R1, R2 = Rvals[-2], Rvals[-1]
        v_rich = (R2*results[R2][amp] - R1*results[R1][amp]) / (R2 - R1)
        print(f"    Richardson 1/R (R={R1},{R2}): {v_rich:+.14f}")

# Fit window sensitivity at R_max=3000
if 3000 in results:
    print(f"\n{'='*78}")
    print("  FIT WINDOW SENSITIVITY (R_max = 3000)")
    print("="*78)

    R_max = 3000
    r_eval = np.linspace(r0, R_max, int(R_max * 20) + 1)
    sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)
    r_arr = sol.t
    u4 = r_arr * sol.y[4]
    u5 = r_arr * sol.y[6]

    windows = [(500,2800), (800,2800), (1000,2800), (1500,2800), (1000,2500)]
    for rmin, rmax in windows:
        B4w, A4w = extract_amplitudes(r_arr, u4, rmin, rmax)
        B5w, A5w = extract_amplitudes(r_arr, u5, rmin, rmax)
        print(f"    [{rmin:4d},{rmax:4d}]: B4={B4w:+.12f}, A4={A4w:+.12f}, B5={B5w:+.12f}, A5={A5w:+.12f}")

# Compute alpha_5
print(f"\n{'='*78}")
print("  ALPHA_5 COMPUTATION")
print("="*78)

if results:
    R_best = max(results.keys())
    B4v = results[R_best]['B4']
    A4v = results[R_best]['A4']
    B5v = results[R_best]['B5']
    A5v = results[R_best]['A5']

    A2v = A2_known
    B2v = B2_known
    B3v = B3_known
    A3v = A3_known

    # alpha_4 = B4 + A2*A3 - B2*A2^2/2
    alpha4 = B4v + A2v*A3v - B2v*A2v**2/2
    print(f"  alpha_4 = {alpha4:+.12f}  (expect ~-0.0246)")

    # alpha_5 = B5 + A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
    t1 = A3v**2/2
    t2 = A2v*A4v
    t3 = -B2v*A2v*A3v
    t4 = A2v**2*(B2v**2-B3v)/2
    t5 = -A2v**4/8
    alpha5 = B5v + t1 + t2 + t3 + t4 + t5

    print(f"\n  alpha_5 terms:")
    print(f"    B5               = {B5v:+.12f}")
    print(f"    A3^2/2           = {t1:+.12f}")
    print(f"    A2*A4            = {t2:+.12f}")
    print(f"    -B2*A2*A3        = {t3:+.12f}")
    print(f"    A2^2*(B2^2-B3)/2 = {t4:+.12f}")
    print(f"    -A2^4/8          = {t5:+.12f}")
    print(f"    --------------------------------")
    print(f"    alpha_5          = {alpha5:+.12f}")
    print(f"    From ODE:        ~ 0.02750")
    print(f"    Difference:      = {alpha5 - 0.02750:+.6f}")

    # eta_sym check
    print(f"\n  eta_sym(d) = 1 + alpha_3*d^2 + alpha_5*d^4 + ...")
    alpha3 = B3v + A2v**2/2
    print(f"    alpha_3 = {alpha3:+.12f}  (known: 0.089722...)")
    print(f"    alpha_5 = {alpha5:+.12f}")

    # Verify: compute eta_sym from perturbation at d=0.1
    d = 0.1
    eta_pert = 1 + alpha3*d**2 + alpha5*d**4
    print(f"\n  eta_sym(0.1) from perturbation: {eta_pert:.12f}")

print(f"\n{'='*78}")
