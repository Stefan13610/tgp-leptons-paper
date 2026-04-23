#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c32_B5_greens_integral.py:
Compute B4, A4, B5, A5 via Green's function integrals.

For f_n'' + (2/r)f_n' + f_n = S_n, the solution u_n = r*f_n satisfies:
  u_n'' + u_n = r*S_n(r)

Green's function gives:
  u_n(r) = int_0^r sin(r-t) * t*S_n(t) dt

As r -> inf, using sin(r-t) = sin(r)cos(t) - cos(r)sin(t):
  B_n = int_0^inf cos(t) * t * S_n(t) dt
  A_n = -int_0^inf sin(t) * t * S_n(t) dt

Strategy:
- Integrate ODE for (f2, f3, f4) with dense output up to R=5000
- Accumulate B_n, A_n integrals using Simpson's rule on fine grid
- Use truncation at different R values for Richardson extrapolation
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

print("="*78)
print("  B4, A4, B5, A5 via Green's function integrals")
print("="*78)

def f1(r):
    if r < 1e-10: return 1.0 - r**2/6 + r**4/120
    return np.sin(r)/r

def f1p(r):
    if r < 1e-10: return -r/3 + r**3/30 - r**5/840
    return np.cos(r)/r - np.sin(r)/r**2

# ODE for f2, f3, f4 (we need these to construct S4, S5)
def rhs(r, y):
    f2, f2p, f3, f3p, f4, f4p = y
    if r < 1e-10:
        return [f2p, 0, f3p, 0, f4p, 0]
    fv = f1(r)
    fp = f1p(r)
    tr = 2.0/r
    s2 = -fp**2
    s3 = fv*fp**2 - 2*fp*f2p
    s4 = 2*fv*fp*f2p + f2*fp**2 - fv**2*fp**2 - 2*fp*f3p - f2p**2
    return [f2p, s2 - tr*f2p - f2,
            f3p, s3 - tr*f3p - f3,
            f4p, s4 - tr*f4p - f4]

r0 = 1e-6
y0 = [0.0]*6

# Step 1: Integrate ODE with very fine output
R_max = 5000
dr = 0.05  # fine enough for Simpson
N = int(R_max / dr)
r_eval = np.linspace(r0, R_max, N + 1)

print(f"\n  Integrating f2, f3, f4 to R_max={R_max} with {N+1} points...")
sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)
print(f"  Integration: {sol.message}")

r_arr = sol.t
f2_arr = sol.y[0]
f2p_arr = sol.y[1]
f3_arr = sol.y[2]
f3p_arr = sol.y[3]
f4_arr = sol.y[4]
f4p_arr = sol.y[5]

# Precompute f1 values
f1_arr = np.array([f1(r) for r in r_arr])
f1p_arr = np.array([f1p(r) for r in r_arr])

# Source terms
S2_arr = -f1p_arr**2
S3_arr = f1_arr*f1p_arr**2 - 2*f1p_arr*f2p_arr
S4_arr = (2*f1_arr*f1p_arr*f2p_arr + f2_arr*f1p_arr**2
          - f1_arr**2*f1p_arr**2 - 2*f1p_arr*f3p_arr - f2p_arr**2)
S5_arr = (-2*f1p_arr*f4p_arr - 2*f2p_arr*f3p_arr
          + 2*f1_arr*f1p_arr*f3p_arr + f1_arr*f2p_arr**2
          + 2*f2_arr*f1p_arr*f2p_arr + f3_arr*f1p_arr**2
          - 2*f1_arr**2*f1p_arr*f2p_arr - 2*f1_arr*f2_arr*f1p_arr**2
          + f1_arr**3*f1p_arr**2)

# Step 2: Compute B_n, A_n integrals via cumulative Simpson
# B_n(R) = int_0^R cos(t)*t*S_n(t) dt
# A_n(R) = -int_0^R sin(t)*t*S_n(t) dt

cos_arr = np.cos(r_arr)
sin_arr = np.sin(r_arr)

def cumulative_simpson(r, y):
    """Cumulative integral using Simpson's rule on uniform-ish grid."""
    n = len(r)
    result = np.zeros(n)
    # Use composite Simpson on pairs of intervals
    for i in range(2, n, 2):
        h = (r[i] - r[i-2]) / 2
        result[i] = result[i-2] + h/3 * (y[i-2] + 4*y[i-1] + y[i])
    # Fill odd indices by interpolation
    for i in range(1, n, 2):
        if i+1 < n:
            result[i] = (result[i-1] + result[i+1]) / 2  # rough, ok for monitoring
        else:
            result[i] = result[i-1]
    return result

# Actually, let's use numpy trapz for different truncation points
def integral_at_R(integrand, R_trunc):
    """Integrate from 0 to R_trunc using trapezoid rule."""
    mask = r_arr <= R_trunc + dr/2
    return np.trapezoid(integrand[mask], r_arr[mask])

# Compute integrands
integrand_B2 = cos_arr * r_arr * S2_arr
integrand_A2 = -sin_arr * r_arr * S2_arr
integrand_B3 = cos_arr * r_arr * S3_arr
integrand_A3 = -sin_arr * r_arr * S3_arr
integrand_B4 = cos_arr * r_arr * S4_arr
integrand_A4 = -sin_arr * r_arr * S4_arr
integrand_B5 = cos_arr * r_arr * S5_arr
integrand_A5 = -sin_arr * r_arr * S5_arr

# Evaluate at multiple truncation points
import math
ln3 = math.log(3)
B2_known = 0.5 - ln3/8
A2_known = math.pi/8
B3_known = 0.012615939290114712
A3_known = -0.215712018451450166

R_values = [500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]

print(f"\n  --- Integrals at different truncation R ---")
print(f"  {'R':>6s}  {'B2':>16s}  {'A2':>16s}  {'B3':>16s}  {'A3':>16s}  {'B4':>16s}  {'A4':>16s}  {'B5':>16s}  {'A5':>16s}")

results = {}
for R in R_values:
    B2 = integral_at_R(integrand_B2, R)
    A2 = integral_at_R(integrand_A2, R)
    B3 = integral_at_R(integrand_B3, R)
    A3 = integral_at_R(integrand_A3, R)
    B4 = integral_at_R(integrand_B4, R)
    A4 = integral_at_R(integrand_A4, R)
    B5 = integral_at_R(integrand_B5, R)
    A5 = integral_at_R(integrand_A5, R)
    results[R] = {'B2':B2,'A2':A2,'B3':B3,'A3':A3,'B4':B4,'A4':A4,'B5':B5,'A5':A5}
    print(f"  {R:6d}  {B2:+16.12f}  {A2:+16.12f}  {B3:+16.12f}  {A3:+16.12f}  "
          f"{B4:+16.12f}  {A4:+16.12f}  {B5:+16.12f}  {A5:+16.12f}")

# Cross-check known values
print(f"\n  --- Cross-check at R=5000 ---")
R = 5000
print(f"    B2 = {results[R]['B2']:+.14f}  (known: {B2_known:+.14f}, err: {results[R]['B2']-B2_known:+.2e})")
print(f"    A2 = {results[R]['A2']:+.14f}  (known: {A2_known:+.14f}, err: {results[R]['A2']-A2_known:+.2e})")
print(f"    B3 = {results[R]['B3']:+.14f}  (known: {B3_known:+.14f}, err: {results[R]['B3']-B3_known:+.2e})")
print(f"    A3 = {results[R]['A3']:+.14f}  (known: {A3_known:+.14f}, err: {results[R]['A3']-A3_known:+.2e})")

# Richardson for B4, A4, B5, A5
print(f"\n{'='*78}")
print("  RICHARDSON EXTRAPOLATION")
print("="*78)

# For conditionally convergent oscillatory integrals, the truncation error
# goes as ~cos(R)/R or ~sin(R)/R. After averaging over many oscillations,
# the running average converges as ~1/R.
# Try Cesaro-style averaging or Richardson

# Simple approach: average consecutive truncations to cancel oscillations
def cesaro_average(R_vals, values):
    """Average over pairs to reduce oscillatory truncation error."""
    avgs = []
    Rs = []
    for i in range(len(R_vals)-1):
        avgs.append((values[i] + values[i+1]) / 2)
        Rs.append((R_vals[i] + R_vals[i+1]) / 2)
    return Rs, avgs

# Better: use truncation at half-integer multiples of pi to cancel cos/sin terms
# Or: fit I(R) = I_inf + a*cos(R)/R + b*sin(R)/R + c/R^2

print(f"\n  Fitting I(R) = I_inf + a*cos(R)/R + b*sin(R)/R + c/R^2:")
from numpy.linalg import lstsq

for amp_name in ['B2', 'A2', 'B3', 'A3', 'B4', 'A4', 'B5', 'A5']:
    R_fit = np.array([R for R in R_values if R >= 1000], dtype=float)
    vals = np.array([results[int(R)][amp_name] for R in R_fit])

    # Basis: 1, cos(R)/R, sin(R)/R, 1/R^2
    M = np.column_stack([
        np.ones_like(R_fit),
        np.cos(R_fit)/R_fit,
        np.sin(R_fit)/R_fit,
        1.0/R_fit**2
    ])

    coeffs, _, _, _ = lstsq(M, vals, rcond=None)
    I_inf = coeffs[0]

    known = {'B2': B2_known, 'A2': A2_known, 'B3': B3_known, 'A3': A3_known}
    if amp_name in known:
        err = I_inf - known[amp_name]
        print(f"    {amp_name}: I_inf = {I_inf:+.14f}  (err: {err:+.2e})")
    else:
        print(f"    {amp_name}: I_inf = {I_inf:+.14f}")

# Also try with more basis functions
print(f"\n  Extended fit: + cos(R)/R^2, sin(R)/R^2:")
for amp_name in ['B4', 'A4', 'B5', 'A5']:
    R_fit = np.array([R for R in R_values if R >= 500], dtype=float)
    vals = np.array([results[int(R)][amp_name] for R in R_fit])

    M = np.column_stack([
        np.ones_like(R_fit),
        np.cos(R_fit)/R_fit,
        np.sin(R_fit)/R_fit,
        1.0/R_fit**2,
        np.cos(R_fit)/R_fit**2,
        np.sin(R_fit)/R_fit**2,
    ])

    coeffs, _, _, _ = lstsq(M, vals, rcond=None)
    I_inf = coeffs[0]
    print(f"    {amp_name}: I_inf = {I_inf:+.14f}")

# Compute alpha_4, alpha_5
print(f"\n{'='*78}")
print("  ALPHA_4, ALPHA_5 from Green's function integrals")
print("="*78)

# Use the fitted I_inf values
R_fit = np.array([R for R in R_values if R >= 1000], dtype=float)

best = {}
for amp_name in ['B4', 'A4', 'B5', 'A5']:
    vals = np.array([results[int(R)][amp_name] for R in R_fit])
    M = np.column_stack([
        np.ones_like(R_fit),
        np.cos(R_fit)/R_fit,
        np.sin(R_fit)/R_fit,
        1.0/R_fit**2
    ])
    coeffs, _, _, _ = lstsq(M, vals, rcond=None)
    best[amp_name] = coeffs[0]

B4v = best['B4']
A4v = best['A4']
B5v = best['B5']
A5v = best['A5']

A2v = A2_known
B2v = B2_known
B3v = B3_known
A3v = A3_known

alpha4 = B4v + A2v*A3v - B2v*A2v**2/2
alpha5 = B5v + A3v**2/2 + A2v*A4v - B2v*A2v*A3v + A2v**2*(B2v**2-B3v)/2 - A2v**4/8

print(f"  B4 = {B4v:+.14f}")
print(f"  A4 = {A4v:+.14f}")
print(f"  B5 = {B5v:+.14f}")
print(f"  A5 = {A5v:+.14f}")
print(f"\n  alpha_4 = {alpha4:+.14f}  (expect ~-0.0246)")
print(f"  alpha_5 = {alpha5:+.14f}  (expect ~0.0275)")

print(f"\n{'='*78}")
