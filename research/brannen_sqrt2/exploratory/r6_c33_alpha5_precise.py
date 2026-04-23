#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c33_alpha5_precise.py:
High-precision alpha_5 via 5th-order perturbation ODE with R_max up to 5000.

Uses tail fitting approach (which gave ~1e-6 errors in known amplitudes at R=3000).
Extends to R=5000 for better convergence.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  HIGH-PRECISION alpha_5 via 5th-order perturbation ODE")
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
    # 4-term fit: B*sin + A*cos + C*sin/r + D*cos/r
    M = np.column_stack([np.sin(r), np.cos(r), np.sin(r)/r, np.cos(r)/r])
    coeffs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    return coeffs[0], coeffs[1]

def extract_amplitudes_6(r_arr, u_arr, r_fit_min, r_fit_max):
    """6-term fit including 1/r^2 corrections."""
    mask = (r_arr >= r_fit_min) & (r_arr <= r_fit_max)
    r = r_arr[mask]
    u = u_arr[mask]
    M = np.column_stack([np.sin(r), np.cos(r),
                         np.sin(r)/r, np.cos(r)/r,
                         np.sin(r)/r**2, np.cos(r)/r**2])
    coeffs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    return coeffs[0], coeffs[1]

# Solve for R_max = 2000, 3000, 5000
R_max_list = [2000, 3000, 5000]
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

    r_arr = sol.t
    u2 = r_arr * sol.y[0]
    u3 = r_arr * sol.y[2]
    u4 = r_arr * sol.y[4]
    u5 = r_arr * sol.y[6]

    # Multiple fit windows
    windows = [
        (R_max*0.3, R_max*0.95),
        (R_max*0.4, R_max*0.95),
        (R_max*0.5, R_max*0.95),
        (R_max*0.4, R_max*0.85),
    ]

    print(f"  4-term fit:")
    for rmin, rmax in windows:
        B4, A4 = extract_amplitudes(r_arr, u4, rmin, rmax)
        B5, A5 = extract_amplitudes(r_arr, u5, rmin, rmax)
        B3, A3 = extract_amplitudes(r_arr, u3, rmin, rmax)
        print(f"    [{rmin:6.0f},{rmax:6.0f}]: B4={B4:+.12f} A4={A4:+.12f} B5={B5:+.12f} "
              f"| B3 err={B3-B3_known:+.1e}, A3 err={A3-A3_known:+.1e}")

    print(f"  6-term fit:")
    for rmin, rmax in windows:
        B4, A4 = extract_amplitudes_6(r_arr, u4, rmin, rmax)
        B5, A5 = extract_amplitudes_6(r_arr, u5, rmin, rmax)
        B3, A3 = extract_amplitudes_6(r_arr, u3, rmin, rmax)
        print(f"    [{rmin:6.0f},{rmax:6.0f}]: B4={B4:+.12f} A4={A4:+.12f} B5={B5:+.12f} "
              f"| B3 err={B3-B3_known:+.1e}, A3 err={A3-A3_known:+.1e}")

    # Store best estimate (6-term, middle window)
    rmin, rmax = R_max*0.4, R_max*0.95
    B2, A2 = extract_amplitudes_6(r_arr, u2, rmin, rmax)
    B3, A3 = extract_amplitudes_6(r_arr, u3, rmin, rmax)
    B4, A4 = extract_amplitudes_6(r_arr, u4, rmin, rmax)
    B5, A5 = extract_amplitudes_6(r_arr, u5, rmin, rmax)
    results[R_max] = {'B2':B2,'A2':A2,'B3':B3,'A3':A3,'B4':B4,'A4':A4,'B5':B5,'A5':A5}

# Final computation
print(f"\n{'='*78}")
print("  FINAL RESULTS")
print("="*78)

# Print convergence
Rvals = sorted(results.keys())
for amp in ['B4', 'A4', 'B5', 'A5']:
    print(f"\n  {amp} convergence (6-term fit):")
    for R in Rvals:
        print(f"    R={R:5d}: {results[R][amp]:+.14f}")
    if len(Rvals) >= 2:
        R1, R2 = Rvals[-2], Rvals[-1]
        v_rich = (R2*results[R2][amp] - R1*results[R1][amp]) / (R2 - R1)
        print(f"    Richardson 1/R: {v_rich:+.14f}")

# Best estimate: R=5000, 6-term fit
R = 5000
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

print(f"\n  BEST ESTIMATES (R_max=5000, 6-term fit):")
print(f"    B4 = {B4v:+.14f}")
print(f"    A4 = {A4v:+.14f}")
print(f"    B5 = {B5v:+.14f}")
print(f"    A5 = {A5v:+.14f}")

print(f"\n  alpha_4 = B4 + A2*A3 - B2*A2^2/2 = {alpha4:+.14f}")

print(f"\n  alpha_5 terms:")
print(f"    B5               = {B5v:+.14f}")
print(f"    A3^2/2           = {t1:+.14f}")
print(f"    A2*A4            = {t2:+.14f}")
print(f"    -B2*A2*A3        = {t3:+.14f}")
print(f"    A2^2*(B2^2-B3)/2 = {t4:+.14f}")
print(f"    -A2^4/8          = {t5:+.14f}")
print(f"    --------------------------------")
print(f"    alpha_5          = {alpha5:+.14f}")

# Cross-check with ODE measurement
print(f"\n  VERIFICATION:")
print(f"    alpha_4: perturbative = {alpha4:+.10f}, ODE ~ -0.0246")
print(f"    alpha_5: perturbative = {alpha5:+.10f}, ODE ~ 0.02750-0.02752")

# Also compute Richardson for alpha_4, alpha_5 directly
print(f"\n  Richardson extrapolation on alpha_4, alpha_5:")
for R1, R2 in [(2000,3000), (3000,5000), (2000,5000)]:
    if R1 in results and R2 in results:
        a4_1 = results[R1]['B4'] + A2v*A3v - B2v*A2v**2/2
        a4_2 = results[R2]['B4'] + A2v*A3v - B2v*A2v**2/2
        a4_r = (R2*a4_2 - R1*a4_1) / (R2 - R1)

        a5_1 = results[R1]['B5'] + t1 + A2v*results[R1]['A4'] + t3 + t4 + t5
        a5_2 = results[R2]['B5'] + t1 + A2v*results[R2]['A4'] + t3 + t4 + t5
        a5_r = (R2*a5_2 - R1*a5_1) / (R2 - R1)

        print(f"    R={R1},{R2}: alpha_4 = {a4_r:+.12f}, alpha_5 = {a5_r:+.12f}")

print(f"\n{'='*78}")
