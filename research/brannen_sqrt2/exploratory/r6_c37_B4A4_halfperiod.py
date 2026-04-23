#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c37_B4A4_halfperiod.py:
High-precision B4, A4, B5, A5 via half-period summation + Wynn-epsilon.

Strategy:
1. Solve coupled f2/f3/f4/f5 ODE with scipy DOP853 + dense_output
2. For each Green's function integral B_n = ∫ cos(t)*t*S_n dt:
   - Core: accurate quadrature on [0, R_start] using scipy.quad
   - Tail: half-period sums a_k = ∫_{k*pi}^{(k+1)*pi} cos(t)*t*S_n dt
   - Accelerate alternating series with Wynn-epsilon
3. Verify against known B2, A2, B3, A3

New analytical results (from r6_c35):
- ∫ sin(t)*t*f1^2*(f1')^2 dt = 5*pi/128
- ∫ sin(t)*t*f1^3*(f1')^2 dt = pi/40
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp, quad
import math

print("="*78)
print("  B4, A4 via half-period summation + Wynn-epsilon")
print("="*78)

# --- ODE setup ---
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
R_max = 5000

print(f"\n  Solving ODE with dense_output to R_max={R_max}...")
sol = solve_ivp(rhs, [r0, R_max], y0, method='DOP853',
                rtol=2.3e-14, atol=1e-20, dense_output=True, max_step=0.1)
print(f"  {sol.message}, nfev={sol.nfev}")

def get_state(r):
    """Return [f2, f2', f3, f3', f4, f4', f5, f5'] at r."""
    if r < r0:
        return [0.0]*8
    return sol.sol(r)

# --- Source terms ---
def S2(r):
    return -f1p(r)**2

def S3(r):
    st = get_state(r)
    return f1(r)*f1p(r)**2 - 2*f1p(r)*st[1]

def S4(r):
    st = get_state(r)
    f2v, f2pv, f3v, f3pv = st[0], st[1], st[2], st[3]
    fp = f1p(r)
    fv = f1(r)
    return 2*fv*fp*f2pv + f2v*fp**2 - fv**2*fp**2 - 2*fp*f3pv - f2pv**2

def S5(r):
    st = get_state(r)
    f2v, f2pv, f3v, f3pv, f4v, f4pv = st[0], st[1], st[2], st[3], st[4], st[5]
    fp = f1p(r)
    fv = f1(r)
    return (-2*fp*f4pv - 2*f2pv*f3pv
            + 2*fv*fp*f3pv + fv*f2pv**2 + 2*f2v*fp*f2pv + f3v*fp**2
            - 2*fv**2*fp*f2pv - 2*fv*f2v*fp**2
            + fv**3*fp**2)

# --- Wynn-epsilon acceleration ---
def wynn_epsilon(partial_sums):
    """Wynn-epsilon acceleration of a sequence of partial sums."""
    n = len(partial_sums)
    if n < 3:
        return partial_sums[-1]

    # Epsilon table
    eps = [[0.0]*(n+1) for _ in range(n+1)]
    for i in range(n):
        eps[i][1] = partial_sums[i]

    for k in range(2, n+1):
        for i in range(n-k+1):
            diff = eps[i+1][k-1] - eps[i][k-1]
            if abs(diff) < 1e-30:
                eps[i][k] = 1e30
            else:
                eps[i][k] = eps[i+1][k-2] + 1.0/diff

    # Best estimate: last even column entry
    best = partial_sums[-1]
    for k in range(1, n+1):
        if k % 2 == 1:  # odd k = even column (eps_0)
            if abs(eps[0][k]) < 1e20:
                best = eps[0][k]

    return best

# --- Half-period summation ---
def halfperiod_integral(integrand_func, R_start_periods=10, n_halfperiods=300):
    """
    Compute ∫_0^∞ integrand(t) dt via:
    1. Core: ∫_0^{R_start} by quadrature
    2. Tail: half-period sums + Wynn-epsilon
    """
    R_start = R_start_periods * math.pi

    # Core integral
    core, core_err = quad(integrand_func, r0, R_start, limit=500, epsabs=1e-16, epsrel=1e-14)

    # Half-period sums
    partial_sums = []
    running = core
    for k in range(n_halfperiods):
        a = R_start + k * math.pi
        b = a + math.pi
        if b > R_max:
            break
        hp, _ = quad(integrand_func, a, b, limit=50, epsabs=1e-16, epsrel=1e-14)
        running += hp
        partial_sums.append(running)

    if len(partial_sums) < 5:
        return running

    # Wynn-epsilon acceleration
    accel = wynn_epsilon(partial_sums)
    return accel

# --- Compute B_n, A_n ---
print(f"\n{'='*78}")
print("  Computing Green's function integrals")
print("="*78)

ln3 = math.log(3)
B2_known = 0.5 - ln3/8
A2_known = math.pi/8
B3_known = 0.012615939290114712
A3_known = -0.215712018451450166

# B2 = ∫ cos(t)*t*S2(t) dt, A2 = -∫ sin(t)*t*S2(t) dt
print("\n  --- S2: -(f1')^2 ---")
B2 = halfperiod_integral(lambda t: math.cos(t)*t*S2(t))
A2 = halfperiod_integral(lambda t: -math.sin(t)*t*S2(t))
print(f"  B2 = {B2:+.15f} (known: {B2_known:+.15f}, err: {B2-B2_known:+.2e})")
print(f"  A2 = {A2:+.15f} (known: {A2_known:+.15f}, err: {A2-A2_known:+.2e})")

# B3
print("\n  --- S3: f1*(f1')^2 - 2*f1'*f2' ---")
B3 = halfperiod_integral(lambda t: math.cos(t)*t*S3(t))
A3 = halfperiod_integral(lambda t: -math.sin(t)*t*S3(t))
print(f"  B3 = {B3:+.15f} (known: {B3_known:+.15f}, err: {B3-B3_known:+.2e})")
print(f"  A3 = {A3:+.15f} (known: {A3_known:+.15f}, err: {A3-A3_known:+.2e})")

# B4, A4
print("\n  --- S4: 2*f1*f1'*f2' + f2*(f1')^2 - f1^2*(f1')^2 - 2*f1'*f3' - (f2')^2 ---")
B4 = halfperiod_integral(lambda t: math.cos(t)*t*S4(t))
A4 = halfperiod_integral(lambda t: -math.sin(t)*t*S4(t))
print(f"  B4 = {B4:+.15f}")
print(f"  A4 = {A4:+.15f}")

# Individual S4 terms for decomposition
print("\n  S4 decomposition:")
B4_T1 = halfperiod_integral(lambda t: math.cos(t)*t*(-f1(t)**2*f1p(t)**2))
A4_T1 = halfperiod_integral(lambda t: -math.sin(t)*t*(-f1(t)**2*f1p(t)**2))
print(f"  T1 (-f1^2*(f1')^2):  B={B4_T1:+.15f}, A={A4_T1:+.15f} (A known: 5pi/128={5*math.pi/128:+.15f})")

B4_rest = B4 - B4_T1
A4_rest = A4 - A4_T1
print(f"  T2+T3+T4+T5:        B={B4_rest:+.15f}, A={A4_rest:+.15f}")

# B5, A5
print("\n  --- S5 ---")
B5 = halfperiod_integral(lambda t: math.cos(t)*t*S5(t))
A5 = halfperiod_integral(lambda t: -math.sin(t)*t*S5(t))
print(f"  B5 = {B5:+.15f}")
print(f"  A5 = {A5:+.15f}")

# --- alpha_4, alpha_5 ---
print(f"\n{'='*78}")
print("  FINAL: alpha_4, alpha_5")
print("="*78)

alpha4 = B4 + A2_known*A3_known - B2_known*A2_known**2/2
t1 = A3_known**2/2
t2 = A2_known*A4
t3 = -B2_known*A2_known*A3_known
t4 = A2_known**2*(B2_known**2-B3_known)/2
t5 = -A2_known**4/8
alpha5 = B5 + t1 + t2 + t3 + t4 + t5

print(f"  B4 = {B4:+.15f}")
print(f"  A4 = {A4:+.15f}")
print(f"  B5 = {B5:+.15f}")
print(f"  A5 = {A5:+.15f}")
print(f"\n  alpha_4 = {alpha4:+.15f}  (expect ~-0.0246)")
print(f"  alpha_5 = {alpha5:+.15f}  (expect ~0.02751)")

# Summary of all amplitudes
print(f"\n  COMPLETE AMPLITUDE TABLE:")
print(f"    B2 = {B2_known:+.15f} (PROVEN: 1/2-ln3/8)")
print(f"    A2 = {A2_known:+.15f} (PROVEN: pi/8)")
print(f"    B3 = {B3_known:+.15f} (30 digits)")
print(f"    A3 = {A3_known:+.15f} (25 digits)")
print(f"    B4 = {B4:+.15f} (this work)")
print(f"    A4 = {A4:+.15f} (this work)")
print(f"    B5 = {B5:+.15f} (this work)")
print(f"    A5 = {A5:+.15f} (this work)")

print(f"\n{'='*78}")
