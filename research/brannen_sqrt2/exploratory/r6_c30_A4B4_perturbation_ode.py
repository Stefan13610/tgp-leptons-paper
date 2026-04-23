#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c30_A4B4_perturbation_ode.py:
Compute A4, B4 (tail amplitudes of u4 = r*f4) via direct ODE integration
of the coupled perturbation equations f2, f3, f4.

Equations (from g'' + (1/g)(g')^2 + (2/r)g' + g = 1, g = 1+h, h = sum d^n f_n):

  f1'' + (2/r)f1' + f1 = 0                          => f1 = sin(r)/r
  f2'' + (2/r)f2' + f2 = -(f1')^2                   =: S2
  f3'' + (2/r)f3' + f3 = f1*(f1')^2 - 2*f1'*f2'     =: S3
  f4'' + (2/r)f4' + f4 = 2*f1*f1'*f2' + f2*(f1')^2
                        - f1^2*(f1')^2 - 2*f1'*f3'
                        - (f2')^2                    =: S4

ICs: f1(0)=1, f1'(0)=0; fn(0)=0, fn'(0)=0 for n>=2.

Strategy:
- Use scipy DOP853 with float64, R_max up to 2000-3000
- Extract tail amplitudes via least-squares fit of u_n = r*f_n
  to B_n*sin(r) + A_n*cos(r) + C_n*sin(r)/r + D_n*cos(r)/r
- Richardson extrapolation for truncation error
- Cross-check: B2, A2, B3, A3 against known values
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp

print("="*78)
print("  A4, B4 via perturbation ODE (scipy DOP853)")
print("="*78)

# f1 and f1' (analytical)
def f1(r):
    if r < 1e-10:
        return 1.0 - r**2/6 + r**4/120
    return np.sin(r)/r

def f1p(r):
    if r < 1e-10:
        return -r/3 + r**3/30 - r**5/840
    return np.cos(r)/r - np.sin(r)/r**2

# Source terms
def S2(r, f1v, f1pv):
    return -f1pv**2

def S3(r, f1v, f1pv, f2pv):
    return f1v * f1pv**2 - 2*f1pv*f2pv

def S4(r, f1v, f1pv, f2v, f2pv, f3pv):
    return (2*f1v*f1pv*f2pv + f2v*f1pv**2
            - f1v**2 * f1pv**2 - 2*f1pv*f3pv
            - f2pv**2)

# ODE system: y = [f2, f2', f3, f3', f4, f4']
def rhs(r, y):
    f2v, f2pv, f3v, f3pv, f4v, f4pv = y
    f1v = f1(r)
    f1pv = f1p(r)

    if r < 1e-10:
        # At r~0, use L'Hopital for (2/r)*f' terms
        # f_n'' + (2/r)*f_n' + f_n = S_n
        # Near r=0: 3*f_n''(0) + f_n(0) = S_n(0) (from regularity)
        # But we start at r=r0 > 0, so this branch rarely triggers
        return [f2pv, 0, f3pv, 0, f4pv, 0]

    two_over_r = 2.0/r

    s2 = S2(r, f1v, f1pv)
    s3 = S3(r, f1v, f1pv, f2pv)
    s4 = S4(r, f1v, f1pv, f2v, f2pv, f3pv)

    f2pp = s2 - two_over_r*f2pv - f2v
    f3pp = s3 - two_over_r*f3pv - f3v
    f4pp = s4 - two_over_r*f4pv - f4v

    return [f2pv, f2pp, f3pv, f3pp, f4pv, f4pp]

# Start from small r0 to avoid 2/r singularity
r0 = 1e-6
# At r=r0, use Taylor series for ICs
# f2(r) ~ a2*r^2 + ... where 3*2*a2 + 0 = S2(0) = -(f1'(0))^2 = 0
# Actually f1'(0) = 0, so S2(0) = 0, S3(0) = 0, S4(0) = 0
# Need higher-order Taylor: f1' ~ -r/3 + r^3/30, so (f1')^2 ~ r^2/9 near 0
# S2 = -r^2/9 + O(r^4)
# For f2: f2'' + (2/r)f2' + f2 = -r^2/9
# Try f2 ~ a*r^2: 2a + 4a + 0 = -r^2/9... no, let's just start at r0 with zeros
y0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Known values for cross-check
import math
ln3 = math.log(3)
B2_known = 0.5 - ln3/8          # I_cos
A2_known = math.pi/8             # pi/8
B3_known = 0.012615939290114712  # P_cos
A3_known = -0.215712018451450166 # -P_sin

def extract_amplitudes(r_arr, u_arr, r_fit_min, r_fit_max):
    """Extract B, A from u ~ B*sin(r) + A*cos(r) + C*sin(r)/r + D*cos(r)/r"""
    mask = (r_arr >= r_fit_min) & (r_arr <= r_fit_max)
    r = r_arr[mask]
    u = u_arr[mask]

    # Basis: sin(r), cos(r), sin(r)/r, cos(r)/r
    M = np.column_stack([
        np.sin(r),
        np.cos(r),
        np.sin(r)/r,
        np.cos(r)/r
    ])

    coeffs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    return coeffs[0], coeffs[1]  # B, A

# Solve for multiple R_max values
R_max_list = [500, 1000, 1500, 2000, 3000]

results = {}
for R_max in R_max_list:
    print(f"\n  --- R_max = {R_max} ---")

    # Dense output for tail extraction
    r_eval = np.linspace(r0, R_max, int(R_max * 20) + 1)

    sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval,
                    max_step=0.1)

    if not sol.success:
        print(f"    Integration failed: {sol.message}")
        continue

    r_arr = sol.t
    f2_arr = sol.y[0]
    f3_arr = sol.y[2]
    f4_arr = sol.y[4]

    # Compute u_n = r * f_n
    u2 = r_arr * f2_arr
    u3 = r_arr * f3_arr
    u4 = r_arr * f4_arr

    # Also compute f1 for u1
    f1_arr = np.array([f1(r) for r in r_arr])
    u1 = r_arr * f1_arr  # should be sin(r)

    # Fit windows
    r_fit_min = R_max * 0.4
    r_fit_max = R_max * 0.95

    B1, A1 = extract_amplitudes(r_arr, u1, r_fit_min, r_fit_max)
    B2, A2 = extract_amplitudes(r_arr, u2, r_fit_min, r_fit_max)
    B3, A3 = extract_amplitudes(r_arr, u3, r_fit_min, r_fit_max)
    B4, A4 = extract_amplitudes(r_arr, u4, r_fit_min, r_fit_max)

    print(f"    B1 = {B1:+.15f}  (expect 1.0)")
    print(f"    A1 = {A1:+.15f}  (expect 0.0)")
    print(f"    B2 = {B2:+.15f}  (known: {B2_known:+.15f}, err: {B2-B2_known:+.2e})")
    print(f"    A2 = {A2:+.15f}  (known: {A2_known:+.15f}, err: {A2-A2_known:+.2e})")
    print(f"    B3 = {B3:+.15f}  (known: {B3_known:+.15f}, err: {B3-B3_known:+.2e})")
    print(f"    A3 = {A3:+.15f}  (known: {A3_known:+.15f}, err: {A3-A3_known:+.2e})")
    print(f"    B4 = {B4:+.15f}")
    print(f"    A4 = {A4:+.15f}")

    results[R_max] = {'B2': B2, 'A2': A2, 'B3': B3, 'A3': A3, 'B4': B4, 'A4': A4}

# Richardson extrapolation
print(f"\n{'='*78}")
print("  RICHARDSON EXTRAPOLATION for A4, B4")
print("="*78)

if len(results) >= 3:
    # Use 1/R^2 Richardson between pairs
    Rvals = sorted(results.keys())

    print(f"\n  Raw B4 values:")
    for R in Rvals:
        print(f"    R={R:5d}: B4 = {results[R]['B4']:+.15f}")

    print(f"\n  Raw A4 values:")
    for R in Rvals:
        print(f"    R={R:5d}: A4 = {results[R]['A4']:+.15f}")

    # 2-point Richardson: f(R) ~ f_inf + C/R
    # f_inf = (R2*f(R2) - R1*f(R1)) / (R2 - R1)
    print(f"\n  Richardson 1/R (2-point):")
    for i in range(len(Rvals)-1):
        R1, R2 = Rvals[i], Rvals[i+1]
        B4_rich = (R2*results[R2]['B4'] - R1*results[R1]['B4']) / (R2 - R1)
        A4_rich = (R2*results[R2]['A4'] - R1*results[R1]['A4']) / (R2 - R1)
        print(f"    R={R1},{R2}: B4_inf = {B4_rich:+.15f}, A4_inf = {A4_rich:+.15f}")

    # 2-point Richardson: f(R) ~ f_inf + C/R^2
    # f_inf = (R2^2*f(R2) - R1^2*f(R1)) / (R2^2 - R1^2)
    print(f"\n  Richardson 1/R^2 (2-point):")
    for i in range(len(Rvals)-1):
        R1, R2 = Rvals[i], Rvals[i+1]
        B4_rich = (R2**2*results[R2]['B4'] - R1**2*results[R1]['B4']) / (R2**2 - R1**2)
        A4_rich = (R2**2*results[R2]['A4'] - R1**2*results[R1]['A4']) / (R2**2 - R1**2)
        print(f"    R={R1},{R2}: B4_inf = {B4_rich:+.15f}, A4_inf = {A4_rich:+.15f}")

    # Also do for B2, A2, B3, A3 as cross-check
    print(f"\n  Cross-check: Richardson 1/R for B2, A2, B3, A3:")
    R1, R2 = Rvals[-2], Rvals[-1]
    for amp_name in ['B2', 'A2', 'B3', 'A3']:
        v1 = results[R1][amp_name]
        v2 = results[R2][amp_name]
        v_rich = (R2*v2 - R1*v1) / (R2 - R1)
        known = {'B2': B2_known, 'A2': A2_known, 'B3': B3_known, 'A3': A3_known}[amp_name]
        print(f"    {amp_name}: rich = {v_rich:+.15f}, known = {known:+.15f}, err = {v_rich-known:+.2e}")

# Also try different fit windows for the largest R_max
if 3000 in results:
    print(f"\n{'='*78}")
    print("  FIT WINDOW SENSITIVITY (R_max = 3000)")
    print("="*78)

    R_max = 3000
    # Re-solve
    r_eval = np.linspace(r0, R_max, int(R_max * 20) + 1)
    sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval,
                    max_step=0.1)
    r_arr = sol.t
    u4 = r_arr * sol.y[4]
    u3 = r_arr * sol.y[2]
    u2 = r_arr * sol.y[0]

    windows = [(500, 2800), (800, 2800), (1000, 2800), (1500, 2800), (1000, 2500)]
    for r_min, r_max_fit in windows:
        B4w, A4w = extract_amplitudes(r_arr, u4, r_min, r_max_fit)
        B3w, A3w = extract_amplitudes(r_arr, u3, r_min, r_max_fit)
        print(f"    [{r_min:4d},{r_max_fit:4d}]: B4={B4w:+.12f}, A4={A4w:+.12f}  "
              f"| B3={B3w:+.12f} (err={B3w-B3_known:+.1e}), A3={A3w:+.12f} (err={A3w-A3_known:+.1e})")

# Compute partial alpha_5 with the measured A4
print(f"\n{'='*78}")
print("  PARTIAL alpha_5 COMPUTATION")
print("="*78)

if results:
    R_best = max(results.keys())
    A4_best = results[R_best]['A4']
    B4_best = results[R_best]['B4']

    A2v = A2_known
    B2v = B2_known
    B3v = B3_known
    A3v = A3_known

    # alpha_5 = B5 + A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
    # Without B5:
    a5_no_B5 = A3v**2/2 + A2v*A4_best - B2v*A2v*A3v + A2v**2*(B2v**2-B3v)/2 - A2v**4/8

    print(f"  Using A4 = {A4_best:+.12f} from R_max = {R_best}")
    print(f"  Terms:")
    print(f"    A3^2/2           = {A3v**2/2:+.12f}")
    print(f"    A2*A4            = {A2v*A4_best:+.12f}")
    print(f"    -B2*A2*A3        = {-B2v*A2v*A3v:+.12f}")
    print(f"    A2^2*(B2^2-B3)/2 = {A2v**2*(B2v**2-B3v)/2:+.12f}")
    print(f"    -A2^4/8          = {-A2v**4/8:+.12f}")
    print(f"    Sum (excl B5)    = {a5_no_B5:+.12f}")
    print(f"  alpha_5 = B5 + {a5_no_B5:+.12f}")
    print(f"  From ODE: alpha_5 ~ 0.02750, so B5 ~ {0.02750 - a5_no_B5:+.6f}")

# alpha_4 computation
print(f"\n{'='*78}")
print("  alpha_4 COMPUTATION")
print("="*78)

if results:
    R_best = max(results.keys())
    B4v = results[R_best]['B4']
    A4v = results[R_best]['A4']

    # alpha_4 = B4 + A2*A3 - B2*A2^2/2
    # (from eta_sym coefficient at d^3 in the odd expansion,
    #  but alpha_4 is the d^3 coeff of eta(d) = 1 + alpha_2*d + alpha_3*d^2 + alpha_4*d^3 + ...)
    # Actually: alpha_4 = B4 + A2*A3 - B2*A2^2/2
    alpha4 = B4v + A2_known*A3_known - B2_known*A2_known**2/2

    print(f"  B4 = {B4v:+.12f}")
    print(f"  A2*A3 = {A2_known*A3_known:+.12f}")
    print(f"  -B2*A2^2/2 = {-B2_known*A2_known**2/2:+.12f}")
    print(f"  alpha_4 = B4 + A2*A3 - B2*A2^2/2 = {alpha4:+.12f}")
    print(f"  Expected from ODE: ~ -0.0246")

print(f"\n{'='*78}")
