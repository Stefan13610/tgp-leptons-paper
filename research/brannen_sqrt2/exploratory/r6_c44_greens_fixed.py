#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c44_greens_fixed.py:
High-precision B4, A4, B5, A5 via on-the-fly Green's function integrals.

FIX from r6_c43: half-period boundaries drifted because steps_per_hp=628
doesn't exactly equal π/dr=628.318... After 953 half-periods, the drift
exceeded 1.5 radians, destroying the alternating sign pattern for Wynn-epsilon.

Fix: detect half-period crossings via floor(r/π) instead of counting steps.
Use composite trapezoidal rule (not Simpson) per half-period — simpler and
more robust for variable-length intervals.

Alternative approach: compute the TOTAL integral via composite Simpson on
the full [0, R_max] grid, then use Wynn-epsilon on properly-aligned
half-period partial sums.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import time
from mpmath import mp, mpf, sin, cos, pi, log, nstr, pslq, zeta, catalan, euler, floor

mp.dps = 30

print("="*78)
print("  FIXED GREEN'S INTEGRALS (dps=30, R=3000, dr=0.005)")
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
def rhs_ode(r, y):
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

print(f"\n  R_max = {R_max}, dr = {float(dr)}, N_steps = {N_steps}")
print(f"  dps = {mp.dps}")

r = mpf('1e-6')
y = [mpf(0)]*8

# Running trapezoidal sum for current half-period
# 8 integrands: B2,A2,B3,A3,B4,A4,B5,A5
hp_trap = [mpf(0)] * 8   # trapezoidal accumulator for current HP
prev_integrands = [mpf(0)] * 8  # previous step integrands

# Running totals
running = [mpf(0)] * 8
partial_sums = [[] for _ in range(8)]

current_hp = 0  # floor(r / pi)

t0 = time.time()

for step in range(N_steps + 1):  # +1 to include the final point
    # Compute integrands
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

    cr = cos(r)
    sr = sin(r)

    integrands = [
        cr * r * s2, -sr * r * s2,
        cr * r * s3, -sr * r * s3,
        cr * r * s4, -sr * r * s4,
        cr * r * s5, -sr * r * s5,
    ]

    # Detect half-period crossing
    new_hp = int(floor(r / pi))
    if new_hp > current_hp and step > 0:
        # We crossed a half-period boundary
        # Finalize current half-period with trapezoidal rule
        for k in range(8):
            running[k] += hp_trap[k]
            partial_sums[k].append(running[k])

        current_hp = new_hp
        # Reset trapezoidal accumulator
        hp_trap = [mpf(0)] * 8

    # Trapezoidal rule: add dr/2 * (f_prev + f_curr) for each step
    if step > 0:
        for k in range(8):
            hp_trap[k] += dr / 2 * (prev_integrands[k] + integrands[k])

    prev_integrands = list(integrands)

    # Advance ODE (if not last step)
    if step < N_steps:
        y = rk4_step(rhs_ode, r, y, dr)
        r += dr

    if step > 0 and step % (N_steps // 20) == 0:
        elapsed = time.time() - t0
        pct = 100 * step / N_steps
        n_hp = len(partial_sums[0])
        print(f"    {pct:.0f}% done (r={float(r):.0f}, {n_hp} hp), elapsed={elapsed:.1f}s")

# Finalize last partial half-period
for k in range(8):
    running[k] += hp_trap[k]

t_total = time.time() - t0
n_hp = len(partial_sums[0])
print(f"\n  Integration complete in {t_total:.1f}s, {n_hp} half-periods")

# --- Wynn-epsilon ---
print(f"\n{'='*78}")
print("  WYNN-EPSILON ACCELERATION")
print("="*78)

labels = ['B2', 'A2', 'B3', 'A3', 'B4', 'A4', 'B5', 'A5']
results_raw = list(running)
results_accel = []

for k in range(8):
    accel = wynn_epsilon(partial_sums[k])
    results_accel.append(accel)
    print(f"  {labels[k]}: raw = {nstr(results_raw[k], 16):>20s}  accel = {nstr(accel, 16):>20s}")

B2, A2, B3, A3, B4, A4, B5, A5 = results_accel

# --- Convergence check ---
print(f"\n  Convergence of Wynn-epsilon:")
for checkpoint in [50, 100, 200, 400, 600, 800, 900]:
    if checkpoint > n_hp: break
    vals = [wynn_epsilon(partial_sums[k][:checkpoint]) for k in range(8)]
    print(f"    {checkpoint:4d} hp: B2={nstr(vals[0],15):>19s}  A2={nstr(vals[1],15):>19s}  B3={nstr(vals[2],15):>19s}  A3={nstr(vals[3],15):>19s}")
    print(f"           B4={nstr(vals[4],15):>19s}  A4={nstr(vals[5],15):>19s}  B5={nstr(vals[6],15):>19s}  A5={nstr(vals[7],15):>19s}")

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
    if abs(err) > 0 and abs(known) > 0:
        n_digits = max(0, -int(float(log(abs(err/known), 10))))
    else:
        n_digits = 99
    print(f"  {name} = {nstr(val, 18)}")
    print(f"  {name}_known = {nstr(known, 18)}")
    print(f"  error = {nstr(err, 5)} (~{n_digits} digits)\n")

# --- alpha_4, alpha_5 ---
print(f"{'='*78}")
print("  ALPHA_4, ALPHA_5")
print("="*78)

A2k = A2_known
B2k = B2_known

alpha4 = B4 + A2k*A3_known - B2k*A2k**2/2
alpha5 = B5 + A3_known**2/2 + A2k*A4 - B2k*A2k*A3_known + A2k**2*(B2k**2-B3_known)/2 - A2k**4/8

print(f"  B4 = {nstr(B4, 18)}")
print(f"  A4 = {nstr(A4, 18)}")
print(f"  B5 = {nstr(B5, 18)}")
print(f"  A5 = {nstr(A5, 18)}")
print(f"\n  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")

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
        ("ln2,ln3", [value, log(2), log(3), mpf(1)]),
        ("pi,ln2,ln3", [value, pi, log(2), log(3), mpf(1)]),
        ("zeta3,pi^2", [value, zeta(3), pi**2, pi, mpf(1)]),
        ("full7", [value, pi, pi**2, pi**3, log(2), log(3), zeta(3), catalan, mpf(1)]),
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
