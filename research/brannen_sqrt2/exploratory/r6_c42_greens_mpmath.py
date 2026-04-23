#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c42_greens_mpmath.py:
High-precision B4, A4, B5, A5 via Green's function integrals
computed with mpmath ODE solution.

Strategy:
1. Solve 8-component ODE with mpmath RK4 at dps=30, dr=0.005, R_max=600
2. Store ALL solution points (120000 × 8 mpf values)
3. Compute Green's function integrals:
   B_n = ∫₀^R cos(t)·t·S_n(t) dt  (via half-period Simpson + Wynn-epsilon)
   A_n = -∫₀^R sin(t)·t·S_n(t) dt
4. Half-period sums naturally handle the oscillatory tail

Key advantage over tail fitting: no O(1/r²) fitting bias.
RK4 at dr=0.005 gives global error O(dr⁴) = O(6e-11) → 10+ digits.
Half-period sums + Wynn-epsilon → fast convergence.

Estimated runtime: ~2 min ODE + ~2 min integrals = ~4 min total.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import time
from mpmath import mp, mpf, sin, cos, pi, log, nstr, pslq, zeta, catalan, euler

mp.dps = 30

print("="*78)
print("  GREEN'S FUNCTION INTEGRALS with mpmath ODE (dps=30)")
print("="*78)

# --- f1, f1' ---
def f1(r):
    if r < mpf('1e-10'):
        return mpf(1) - r**2/6 + r**4/120
    return sin(r)/r

def f1p(r):
    if r < mpf('1e-10'):
        return -r/3 + r**3/30 - r**5/840
    return cos(r)/r - sin(r)/r**2

# --- RK4 ---
def rk4_step(rhs_func, r, y, h):
    k1 = rhs_func(r, y)
    k2 = rhs_func(r + h/2, [y[i] + h/2*k1[i] for i in range(len(y))])
    k3 = rhs_func(r + h/2, [y[i] + h/2*k2[i] for i in range(len(y))])
    k4 = rhs_func(r + h, [y[i] + h*k3[i] for i in range(len(y))])
    return [y[i] + h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(len(y))]

# --- RHS ---
def rhs(r, y):
    f2, f2p, f3, f3p, f4, f4p, f5, f5p = y
    fv = f1(r)
    fp = f1p(r)
    if r < mpf('1e-10'):
        return [f2p, mpf(0), f3p, mpf(0), f4p, mpf(0), f5p, mpf(0)]
    tr = mpf(2)/r
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

# --- ODE Integration, store all ---
R_max = 600
dr = mpf('0.005')
N_steps = int(R_max / float(dr))

print(f"\n  R_max = {R_max}, dr = {float(dr)}, N_steps = {N_steps}")
print(f"  dps = {mp.dps}")

r_grid = []
y_grid = []  # store [f2, f2', f3, f3', f4, f4', f5, f5']

r = mpf('1e-6')
y = [mpf(0)]*8

t0 = time.time()

# Store initial point
r_grid.append(r)
y_grid.append(list(y))

for step in range(N_steps):
    y = rk4_step(rhs, r, y, dr)
    r += dr
    r_grid.append(r)
    y_grid.append(list(y))
    if (step + 1) % (N_steps // 10) == 0:
        elapsed = time.time() - t0
        pct = 100*(step+1)/N_steps
        print(f"    {pct:.0f}% done (r={float(r):.0f}), elapsed={elapsed:.1f}s")

t_total = time.time() - t0
print(f"\n  ODE integration complete in {t_total:.1f}s, {len(r_grid)} points stored")

# --- Source terms at grid points ---
def S2_at(i):
    rv = r_grid[i]
    return -f1p(rv)**2

def S3_at(i):
    rv = r_grid[i]
    f2p_val = y_grid[i][1]
    return f1(rv)*f1p(rv)**2 - 2*f1p(rv)*f2p_val

def S4_at(i):
    rv = r_grid[i]
    f2_val, f2p_val = y_grid[i][0], y_grid[i][1]
    f3p_val = y_grid[i][3]
    fv, fp = f1(rv), f1p(rv)
    return 2*fv*fp*f2p_val + f2_val*fp**2 - fv**2*fp**2 - 2*fp*f3p_val - f2p_val**2

def S5_at(i):
    rv = r_grid[i]
    f2_val, f2p_val = y_grid[i][0], y_grid[i][1]
    f3_val, f3p_val = y_grid[i][2], y_grid[i][3]
    f4p_val = y_grid[i][5]
    fv, fp = f1(rv), f1p(rv)
    return (-2*fp*f4p_val - 2*f2p_val*f3p_val
            + 2*fv*fp*f3p_val + fv*f2p_val**2 + 2*f2_val*fp*f2p_val + f3_val*fp**2
            - 2*fv**2*fp*f2p_val - 2*fv*f2_val*fp**2
            + fv**3*fp**2)

# --- Composite Simpson on each half-period ---
def simpson_halfperiod(integrand_B, integrand_A, start_idx, end_idx):
    """Simpson's rule for ∫ integrand(t) dt on [r_grid[start], r_grid[end]].
    Uses the stored grid points. Requires even number of subintervals."""
    n = end_idx - start_idx
    if n <= 0: return mpf(0), mpf(0)
    if n % 2 == 1:
        end_idx -= 1
        n -= 1
    if n == 0: return mpf(0), mpf(0)

    h = dr
    sum_B = integrand_B(start_idx) + integrand_B(end_idx)
    sum_A = integrand_A(start_idx) + integrand_A(end_idx)
    for j in range(1, n, 2):
        sum_B += 4 * integrand_B(start_idx + j)
        sum_A += 4 * integrand_A(start_idx + j)
    for j in range(2, n, 2):
        sum_B += 2 * integrand_B(start_idx + j)
        sum_A += 2 * integrand_A(start_idx + j)

    return sum_B * h/3, sum_A * h/3

# --- Half-period summation + Wynn-epsilon ---
def wynn_epsilon(partial_sums):
    n = len(partial_sums)
    if n < 3: return partial_sums[-1]
    eps = [[mpf(0)]*(n+1) for _ in range(n+1)]
    for i in range(n):
        eps[i][1] = partial_sums[i]
    for k in range(2, n+1):
        for i in range(n-k+1):
            diff = eps[i+1][k-1] - eps[i][k-1]
            if abs(diff) < mpf('1e-50'):
                eps[i][k] = mpf('1e50')
            else:
                eps[i][k] = eps[i+1][k-2] + 1/diff
    best = partial_sums[-1]
    for k in range(1, n+1):
        if k % 2 == 1 and abs(eps[0][k]) < mpf('1e40'):
            best = eps[0][k]
    return best

def compute_BA(S_func, name):
    """Compute B_n = ∫ cos(t)*t*S_n(t) dt, A_n = -∫ sin(t)*t*S_n(t) dt
    via half-period summation + Wynn-epsilon."""

    print(f"\n  Computing {name}...")

    # Steps per half-period: π/dr ≈ 628
    steps_per_hp = int(round(float(pi / dr)))

    # Start from index 1 (skip r ≈ 0)
    # Core: first few half-periods
    R_start_hp = 5  # start half-period summation after 5 half-periods
    core_end_idx = R_start_hp * steps_per_hp
    if core_end_idx >= len(r_grid):
        core_end_idx = len(r_grid) - 1

    # Core integral [0, R_start_hp*π] via Simpson
    def integrand_B(i):
        rv = r_grid[i]
        return cos(rv) * rv * S_func(i)
    def integrand_A(i):
        rv = r_grid[i]
        return -sin(rv) * rv * S_func(i)

    core_B, core_A = simpson_halfperiod(integrand_B, integrand_A, 0, core_end_idx)

    # Half-period sums
    partial_sums_B = []
    partial_sums_A = []
    running_B = core_B
    running_A = core_A

    max_hp = int(R_max * float(1/pi)) + 1
    for k in range(R_start_hp, max_hp):
        hp_start = k * steps_per_hp
        hp_end = (k + 1) * steps_per_hp
        if hp_end >= len(r_grid):
            hp_end = len(r_grid) - 1
        if hp_start >= len(r_grid) - 1:
            break

        hp_B, hp_A = simpson_halfperiod(integrand_B, integrand_A, hp_start, hp_end)
        running_B += hp_B
        running_A += hp_A
        partial_sums_B.append(running_B)
        partial_sums_A.append(running_A)

    # Wynn-epsilon acceleration
    if len(partial_sums_B) >= 5:
        B_accel = wynn_epsilon(partial_sums_B)
        A_accel = wynn_epsilon(partial_sums_A)
    else:
        B_accel = running_B
        A_accel = running_A

    # Also report raw running sum for comparison
    print(f"    {name}: B_raw = {nstr(running_B, 18)}, B_accel = {nstr(B_accel, 18)}")
    print(f"    {name}: A_raw = {nstr(running_A, 18)}, A_accel = {nstr(A_accel, 18)}")
    print(f"    {len(partial_sums_B)} half-periods used")

    return B_accel, A_accel

# --- Compute all ---
print(f"\n{'='*78}")
print("  GREEN'S FUNCTION INTEGRALS")
print("="*78)

t0 = time.time()

B2, A2 = compute_BA(S2_at, "S2")
B3, A3 = compute_BA(S3_at, "S3")
B4, A4 = compute_BA(S4_at, "S4")
B5, A5 = compute_BA(S5_at, "S5")

t_int = time.time() - t0
print(f"\n  Integral computation: {t_int:.1f}s")

# --- Verification ---
print(f"\n{'='*78}")
print("  VERIFICATION")
print("="*78)

ln3 = log(3)
B2_known = mpf(1)/2 - ln3/8
A2_known = pi/8
B3_known = mpf('0.012615939290114712')
A3_known = mpf('-0.215712018451450166')

print(f"  B2 = {nstr(B2, 20)}")
print(f"  B2_known = {nstr(B2_known, 20)}")
print(f"  B2 error = {nstr(B2 - B2_known, 5)}")
print(f"  A2 = {nstr(A2, 20)}")
print(f"  A2_known = {nstr(A2_known, 20)}")
print(f"  A2 error = {nstr(A2 - A2_known, 5)}")
print(f"  B3 = {nstr(B3, 20)}")
print(f"  B3_known = {nstr(B3_known, 20)}")
print(f"  B3 error = {nstr(B3 - B3_known, 5)}")
print(f"  A3 = {nstr(A3, 20)}")
print(f"  A3_known = {nstr(A3_known, 20)}")
print(f"  A3 error = {nstr(A3 - A3_known, 5)}")

# --- alpha_4, alpha_5 ---
print(f"\n{'='*78}")
print("  ALPHA_4, ALPHA_5")
print("="*78)

# Use known A2, B2, B3, A3 for maximum precision
A2k = A2_known
B2k = B2_known

alpha4 = B4 + A2k*A3_known - B2k*A2k**2/2
alpha5 = B5 + A3_known**2/2 + A2k*A4 - B2k*A2k*A3_known + A2k**2*(B2k**2 - B3_known)/2 - A2k**4/8

print(f"  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")

# Using computed A2, B2 as consistency check
alpha4_raw = B4 + A2*A3 - B2*A2**2/2
alpha5_raw = B5 + A3**2/2 + A2*A4 - B2*A2*A3 + A2**2*(B2**2 - B3)/2 - A2**4/8
print(f"  alpha_4 (raw) = {nstr(alpha4_raw, 18)}")
print(f"  alpha_5 (raw) = {nstr(alpha5_raw, 18)}")

# --- PSLQ ---
print(f"\n{'='*78}")
print("  PSLQ SEARCHES")
print("="*78)

def try_pslq(value, name, max_coeff=1000):
    print(f"\n  --- {name} = {nstr(value, 18)} ---")
    for basis_name, basis in [
        ("pi", [value, pi, mpf(1)]),
        ("pi^2", [value, pi**2, pi, mpf(1)]),
        ("ln", [value, log(2), log(3), mpf(1)]),
        ("pi+ln", [value, pi, log(2), log(3), mpf(1)]),
        ("zeta3", [value, zeta(3), pi**2, pi, mpf(1)]),
        ("full", [value, pi, pi**2, pi**3, log(2), log(3), zeta(3), catalan, mpf(1)]),
        ("cross", [value, pi, pi**2, log(2), log(3), log(2)*log(3), zeta(3), mpf(1)]),
    ]:
        try:
            rel = pslq(basis, maxcoeff=max_coeff)
            if rel:
                print(f"    {basis_name}: {rel}")
        except:
            pass

try_pslq(B4, "B4")
try_pslq(A4, "A4")
try_pslq(B5, "B5")
try_pslq(A5, "A5")
try_pslq(alpha4, "alpha_4")
try_pslq(alpha5, "alpha_5")

# Decomposed quantities
try_pslq(A4 - 5*pi/128, "A4 - 5pi/128")
try_pslq(A5 + pi/40, "A5 + pi/40")

print(f"\n{'='*78}")
print("  FINAL SUMMARY")
print("="*78)
print(f"  B2 = {nstr(B2, 18)} (err: {nstr(B2-B2_known, 3)})")
print(f"  A2 = {nstr(A2, 18)} (err: {nstr(A2-A2_known, 3)})")
print(f"  B3 = {nstr(B3, 18)} (err: {nstr(B3-B3_known, 3)})")
print(f"  A3 = {nstr(A3, 18)} (err: {nstr(A3-A3_known, 3)})")
print(f"  B4 = {nstr(B4, 18)}")
print(f"  A4 = {nstr(A4, 18)}")
print(f"  B5 = {nstr(B5, 18)}")
print(f"  A5 = {nstr(A5, 18)}")
print(f"  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")
print(f"\n{'='*78}")
