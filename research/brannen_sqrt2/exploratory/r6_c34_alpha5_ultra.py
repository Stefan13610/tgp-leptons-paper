#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c34_alpha5_ultra.py:
Ultra-precision alpha_5 via 5th-order perturbation ODE with R_max up to 10000.

At R_max=5000, B3 errors were ~1e-7 with 4-term fit.
R_max=10000 should give ~5e-8, allowing 7+ digits of B4, A4, B5.

Also: decompose S4 into analytical (f1-only) and numerical (f2/f3-dependent)
terms for cross-check.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  ULTRA-PRECISION alpha_5 (R_max up to 10000)")
print("="*78)

def f1(r):
    if r < 1e-10: return 1.0 - r**2/6 + r**4/120
    return np.sin(r)/r

def f1p(r):
    if r < 1e-10: return -r/3 + r**3/30 - r**5/840
    return np.cos(r)/r - np.sin(r)/r**2

def rhs(r, y):
    f2, f2p, f3, f3p, f4, f4p, f5, f5p = y
    if r < 1e-10:
        return [f2p, 0, f3p, 0, f4p, 0, f5p, 0]
    fv = f1(r)
    fp = f1p(r)
    tr = 2.0/r
    s2 = -fp**2
    s3 = fv*fp**2 - 2*fp*f2p
    s4 = 2*fv*fp*f2p + f2*fp**2 - fv**2*fp**2 - 2*fp*f3p - f2p**2
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

# Run for R_max = 5000, 8000, 10000
R_max_list = [5000, 8000, 10000]
results = {}

for R_max in R_max_list:
    print(f"\n  --- R_max = {R_max} ---")
    n_pts = R_max * 20 + 1
    r_eval = np.linspace(r0, R_max, n_pts)

    sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)

    if not sol.success:
        print(f"    Failed: {sol.message}")
        continue
    print(f"    nfev = {sol.nfev}")

    r_arr = sol.t
    u2 = r_arr * sol.y[0]
    u3 = r_arr * sol.y[2]
    u4 = r_arr * sol.y[4]
    u5 = r_arr * sol.y[6]

    # Best fit window: [0.4*R, 0.95*R]
    rmin = R_max * 0.4
    rmax = R_max * 0.95

    B2, A2 = extract_amplitudes(r_arr, u2, rmin, rmax)
    B3, A3 = extract_amplitudes(r_arr, u3, rmin, rmax)
    B4, A4 = extract_amplitudes(r_arr, u4, rmin, rmax)
    B5, A5 = extract_amplitudes(r_arr, u5, rmin, rmax)

    print(f"    B2 err = {B2-B2_known:+.2e}, A2 err = {A2-A2_known:+.2e}")
    print(f"    B3 err = {B3-B3_known:+.2e}, A3 err = {A3-A3_known:+.2e}")
    print(f"    B4 = {B4:+.14f}")
    print(f"    A4 = {A4:+.14f}")
    print(f"    B5 = {B5:+.14f}")
    print(f"    A5 = {A5:+.14f}")

    results[R_max] = {'B2':B2,'A2':A2,'B3':B3,'A3':A3,'B4':B4,'A4':A4,'B5':B5,'A5':A5}

    # Also try multiple windows for stability assessment
    windows = [
        (R_max*0.3, R_max*0.95),
        (R_max*0.4, R_max*0.95),
        (R_max*0.5, R_max*0.95),
        (R_max*0.5, R_max*0.90),
    ]
    print(f"    Window sensitivity:")
    for rm, rx in windows:
        b4, a4 = extract_amplitudes(r_arr, u4, rm, rx)
        b5, a5 = extract_amplitudes(r_arr, u5, rm, rx)
        b3, a3 = extract_amplitudes(r_arr, u3, rm, rx)
        print(f"      [{rm:6.0f},{rx:6.0f}]: B4={b4:+.12f} A4={a4:+.12f} B5={b5:+.12f} "
              f"B3e={b3-B3_known:+.1e}")

# Richardson extrapolation
print(f"\n{'='*78}")
print("  RICHARDSON EXTRAPOLATION")
print("="*78)

Rvals = sorted(results.keys())
for amp in ['B4', 'A4', 'B5', 'A5']:
    vals = [(R, results[R][amp]) for R in Rvals]
    print(f"\n  {amp}:")
    for R, v in vals:
        print(f"    R={R:5d}: {v:+.14f}")
    # 1/R Richardson between last two
    if len(vals) >= 2:
        R1, v1 = vals[-2]
        R2, v2 = vals[-1]
        vr = (R2*v2 - R1*v1) / (R2 - R1)
        print(f"    Richardson 1/R ({R1},{R2}): {vr:+.14f}")
    # 3-point Richardson (if 3 values)
    if len(vals) >= 3:
        R1, v1 = vals[0]
        R2, v2 = vals[1]
        R3, v3 = vals[2]
        # Assuming error ~ C/R: fit v = v_inf + C/R
        # v1 = v_inf + C/R1, v2 = v_inf + C/R2
        # Two unknowns: use least squares with 3 points
        M = np.column_stack([np.ones(3), 1.0/np.array([R1,R2,R3])])
        coeffs, _, _, _ = np.linalg.lstsq(M, np.array([v1,v2,v3]), rcond=None)
        print(f"    3-point fit v_inf: {coeffs[0]:+.14f}")

# Final alpha computation
print(f"\n{'='*78}")
print("  FINAL RESULTS")
print("="*78)

R = max(results.keys())
B4v = results[R]['B4']
A4v = results[R]['A4']
B5v = results[R]['B5']
A5v = results[R]['A5']

A2v = A2_known
B2v = B2_known
B3v = B3_known
A3v = A3_known

alpha4 = B4v + A2v*A3v - B2v*A2v**2/2
t1 = A3v**2/2
t2 = A2v*A4v
t3 = -B2v*A2v*A3v
t4 = A2v**2*(B2v**2-B3v)/2
t5 = -A2v**4/8
alpha5 = B5v + t1 + t2 + t3 + t4 + t5

print(f"  R_max = {R}")
print(f"  B4 = {B4v:+.14f}")
print(f"  A4 = {A4v:+.14f}")
print(f"  B5 = {B5v:+.14f}")
print(f"  A5 = {A5v:+.14f}")
print(f"  alpha_4 = {alpha4:+.14f}")
print(f"  alpha_5 = {alpha5:+.14f}")

# Compare with Richardson-extrapolated values
print(f"\n  Richardson-extrapolated alpha_5:")
for i in range(len(Rvals)-1):
    R1, R2 = Rvals[i], Rvals[i+1]
    a5_1 = results[R1]['B5'] + t1 + A2v*results[R1]['A4'] + t3 + t4 + t5
    a5_2 = results[R2]['B5'] + t1 + A2v*results[R2]['A4'] + t3 + t4 + t5
    a5_r = (R2*a5_2 - R1*a5_1) / (R2 - R1)
    print(f"    R={R1},{R2}: alpha_5 = {a5_r:+.14f}")

print(f"\n{'='*78}")
