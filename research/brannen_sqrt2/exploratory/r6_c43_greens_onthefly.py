#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c43_greens_onthefly.py:
High-precision B4, A4, B5, A5 via on-the-fly Green's function integrals.

Key insight from r6_c42: f1-only integrals (B2, A2) achieved 12-18 digits
with R_max=600, but f2-dependent integrals (A3) had only 4 digits because
the tail beyond R=600 wasn't well captured by Wynn-epsilon.

Solution: integrate ODE and Green's function integrals simultaneously.
- mpmath RK4, dps=30, dr=0.005 (ODE global error ~ 6e-11)
- R_max=3000 (5× more range than r6_c42)
- Compute integrals on-the-fly: no memory overhead
- Half-period accumulation + Wynn-epsilon

Memory: only stores ~1000 half-period partial sums (constant memory!)
Runtime: ~5 min ODE + integrals computed in parallel.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import time
from mpmath import mp, mpf, sin, cos, pi, log, nstr, pslq, zeta, catalan, euler

mp.dps = 30

print("="*78)
print("  ON-THE-FLY GREEN'S INTEGRALS (dps=30, R=3000, dr=0.005)")
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

# --- Source terms from state ---
def get_sources(r, y):
    """Return S2, S3, S4, S5 from r and state y."""
    f2, f2p, f3, f3p, f4, f4p, f5, f5p = y
    fv = f1(r)
    fp = f1p(r)
    s2 = -fp**2
    s3 = fv*fp**2 - 2*fp*f2p
    s4 = 2*fv*fp*f2p + f2*fp**2 - fv**2*fp**2 - 2*fp*f3p - f2p**2
    s5 = (-2*fp*f4p - 2*f2p*f3p
          + 2*fv*fp*f3p + fv*f2p**2 + 2*f2*fp*f2p + f3*fp**2
          - 2*fv**2*fp*f2p - 2*fv*f2*fp**2
          + fv**3*fp**2)
    return s2, s3, s4, s5

# --- Wynn-epsilon ---
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

# --- Integration ---
R_max = 3000
dr = mpf('0.005')
N_steps = int(R_max / float(dr))
steps_per_hp = int(round(float(pi / dr)))  # ~628 steps per half-period

print(f"\n  R_max = {R_max}, dr = {float(dr)}, N_steps = {N_steps}")
print(f"  dps = {mp.dps}")
print(f"  Steps per half-period: {steps_per_hp}")

# --- Simpson accumulators for current half-period ---
# For each source (S2-S5) and each weight (cos*t, -sin*t):
# 8 integrands total: B2,A2,B3,A3,B4,A4,B5,A5

# Accumulate Simpson sums within each half-period
# Simpson: ∫ f dx ≈ h/3 * (f0 + 4f1 + 2f2 + 4f3 + ... + fn)
# Weight pattern: 1, 4, 2, 4, 2, ..., 4, 1

hp_idx = 0  # current half-period index
step_in_hp = 0  # step within current half-period
hp_sum = [mpf(0)] * 8  # running Simpson sum for current half-period

# Running totals and partial sums for Wynn
running = [mpf(0)] * 8
partial_sums = [[] for _ in range(8)]  # list of partial sums for each integral

r = mpf('1e-6')
y = [mpf(0)] * 8

t0 = time.time()

for step in range(N_steps):
    # Compute integrands at current point
    s2, s3, s4, s5 = get_sources(r, y)
    cr = cos(r)
    sr = sin(r)

    integrands = [
        cr * r * s2,      # B2
        -sr * r * s2,     # A2
        cr * r * s3,      # B3
        -sr * r * s3,     # A3
        cr * r * s4,      # B4
        -sr * r * s4,     # A4
        cr * r * s5,      # B5
        -sr * r * s5,     # A5
    ]

    # Simpson weight
    if step_in_hp == 0 or step_in_hp == steps_per_hp:
        w = 1
    elif step_in_hp % 2 == 1:
        w = 4
    else:
        w = 2

    for k in range(8):
        hp_sum[k] += w * integrands[k]

    # Check if we completed a half-period
    if step_in_hp == steps_per_hp:
        # Finalize Simpson for this half-period
        for k in range(8):
            running[k] += hp_sum[k] * dr / 3
            partial_sums[k].append(running[k])

        # Reset for next half-period
        hp_idx += 1
        step_in_hp = 0
        hp_sum = [mpf(0)] * 8
        # Re-add current point as start of new half-period
        for k in range(8):
            hp_sum[k] = integrands[k]  # weight = 1
    else:
        step_in_hp += 1

    # Advance ODE
    y = rk4_step(rhs, r, y, dr)
    r += dr

    if (step + 1) % (N_steps // 20) == 0:
        elapsed = time.time() - t0
        pct = 100 * (step + 1) / N_steps
        n_hp = len(partial_sums[0])
        print(f"    {pct:.0f}% done (r={float(r):.0f}, {n_hp} half-periods), elapsed={elapsed:.1f}s")

t_total = time.time() - t0
n_hp = len(partial_sums[0])
print(f"\n  Integration complete in {t_total:.1f}s, {n_hp} half-periods")

# --- Wynn-epsilon acceleration ---
print(f"\n{'='*78}")
print("  WYNN-EPSILON ACCELERATION")
print("="*78)

labels = ['B2', 'A2', 'B3', 'A3', 'B4', 'A4', 'B5', 'A5']
results_raw = [running[k] for k in range(8)]
results_accel = []

for k in range(8):
    accel = wynn_epsilon(partial_sums[k])
    results_accel.append(accel)
    print(f"  {labels[k]}: raw = {nstr(results_raw[k], 16)}, accel = {nstr(accel, 16)}")

B2, A2, B3, A3, B4, A4, B5, A5 = results_accel

# --- Convergence check: Wynn at intermediate points ---
print(f"\n  Convergence of Wynn-epsilon:")
for checkpoint in [100, 200, 400, 600, 800]:
    if checkpoint > n_hp: break
    vals = [wynn_epsilon(partial_sums[k][:checkpoint]) for k in range(8)]
    print(f"    {checkpoint} hp: B4={nstr(vals[4],14)}, A4={nstr(vals[5],14)}, B5={nstr(vals[6],14)}, A5={nstr(vals[7],14)}")

# --- Verification ---
print(f"\n{'='*78}")
print("  VERIFICATION")
print("="*78)

ln3 = log(3)
B2_known = mpf(1)/2 - ln3/8
A2_known = pi/8
B3_known = mpf('0.012615939290114712')
A3_known = mpf('-0.215712018451450166')

for name, val, known in [
    ('B2', B2, B2_known), ('A2', A2, A2_known),
    ('B3', B3, B3_known), ('A3', A3, A3_known),
]:
    err = val - known
    n_digits = -int(float(log(abs(err/known), 10))) if abs(err) > 0 else 99
    print(f"  {name} = {nstr(val, 18)}")
    print(f"  {name}_known = {nstr(known, 18)}")
    print(f"  error = {nstr(err, 5)} (~{n_digits} digits)\n")

# --- alpha_4, alpha_5 ---
print(f"{'='*78}")
print("  ALPHA_4, ALPHA_5")
print("="*78)

A2k = A2_known
B2k = B2_known

alpha4 = B4 + A2k * A3_known - B2k * A2k**2 / 2
alpha5 = B5 + A3_known**2/2 + A2k*A4 - B2k*A2k*A3_known + A2k**2*(B2k**2 - B3_known)/2 - A2k**4/8

print(f"  B4 = {nstr(B4, 18)}")
print(f"  A4 = {nstr(A4, 18)}")
print(f"  B5 = {nstr(B5, 18)}")
print(f"  A5 = {nstr(A5, 18)}")
print(f"\n  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")

# Raw consistency check
alpha4_raw = B4 + A2*A3 - B2*A2**2/2
alpha5_raw = B5 + A3**2/2 + A2*A4 - B2*A2*A3 + A2**2*(B2**2-B3)/2 - A2**4/8
print(f"  alpha_4 (raw) = {nstr(alpha4_raw, 18)}")
print(f"  alpha_5 (raw) = {nstr(alpha5_raw, 18)}")

# --- PSLQ ---
print(f"\n{'='*78}")
print("  PSLQ SEARCHES")
print("="*78)

def try_pslq(value, name, max_coeff=1000):
    print(f"\n  --- {name} = {nstr(value, 18)} ---")
    bases = [
        ("pi", [value, pi, mpf(1)]),
        ("pi^2", [value, pi**2, pi, mpf(1)]),
        ("pi^2,pi^3", [value, pi**2, pi**3, mpf(1)]),
        ("ln2,ln3", [value, log(2), log(3), mpf(1)]),
        ("pi,ln2,ln3", [value, pi, log(2), log(3), mpf(1)]),
        ("zeta3,pi^2", [value, zeta(3), pi**2, pi, mpf(1)]),
        ("full7", [value, pi, pi**2, log(2), log(3), zeta(3), catalan, mpf(1)]),
        ("pi*ln", [value, pi*log(2), pi*log(3), pi, mpf(1)]),
    ]
    for basis_name, basis in bases:
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
try_pslq(A4 - 5*pi/128, "A4 - 5pi/128")
try_pslq(A5 + pi/40, "A5 + pi/40")

print(f"\n{'='*78}")
print("  FINAL SUMMARY")
print("="*78)
print(f"  B4 = {nstr(B4, 16)}")
print(f"  A4 = {nstr(A4, 16)}")
print(f"  B5 = {nstr(B5, 16)}")
print(f"  A5 = {nstr(A5, 16)}")
print(f"  alpha_4 = {nstr(alpha4, 16)}")
print(f"  alpha_5 = {nstr(alpha5, 16)}")
print(f"\n  Estimated precision: determined by B3/A3 errors")
print(f"  B3 error gave ~{-int(float(log(abs(B3-B3_known)/abs(B3_known), 10)))} digits")
print(f"  A3 error gave ~{-int(float(log(abs(A3-A3_known)/abs(A3_known), 10)))} digits")
print(f"\n{'='*78}")
