#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c41_B4A4_mpmath_longrange.py:
High-precision B4, A4, B5, A5 using mpmath RK4 at dps=25, R_max=3000.

The r6_c40 run showed R_max=500 gives only ~4 digits due to O(1/r) tail
corrections. Increasing to R_max=3000 with dr=0.02 should give 8-10 digits
from the 4-term fit, while keeping runtime under 3 minutes.

Strategy: mpmath RK4 -> 4-term/6-term tail fit -> PSLQ
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import time
from mpmath import mp, mpf, sin, cos, sqrt, pi, log, matrix, lu_solve
from mpmath import nstr, pslq, zeta, catalan, euler as euler_gamma

mp.dps = 25

print("="*78)
print("  HIGH-PRECISION B4/A4/B5/A5 via mpmath RK4 (dps=25, R_max=3000)")
print("="*78)

# --- f1 and f1' (exact) ---
def f1(r):
    if r < mpf('1e-10'):
        return mpf(1) - r**2/6 + r**4/120
    return sin(r)/r

def f1p(r):
    if r < mpf('1e-10'):
        return -r/3 + r**3/30 - r**5/840
    return cos(r)/r - sin(r)/r**2

# --- RK4 step ---
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

# --- Integrate ---
R_max = 3000
dr = mpf('0.02')
N_steps = int(R_max / float(dr))

print(f"\n  R_max = {R_max}, dr = {float(dr)}, N_steps = {N_steps}")
print(f"  dps = {mp.dps}")

r = mpf('1e-6')
y = [mpf(0)]*8

t0 = time.time()

# Save at every integer r
r_save = []
y_save = []
save_interval = int(1.0 / float(dr))  # 50 steps per unit

for step in range(N_steps):
    y = rk4_step(rhs, r, y, dr)
    r += dr
    if (step + 1) % save_interval == 0:
        r_save.append(r)
        y_save.append(list(y))
    if (step + 1) % (N_steps // 10) == 0:
        elapsed = time.time() - t0
        pct = 100*(step+1)/N_steps
        print(f"    {pct:.0f}% done (r={float(r):.0f}), elapsed={elapsed:.1f}s")

t_total = time.time() - t0
print(f"\n  ODE integration complete in {t_total:.1f}s")
print(f"  Saved {len(r_save)} points")

# --- Tail amplitude extraction ---
def extract_4term(r_arr, y_arr, ci, fs, fe):
    """4-term fit: u = B*sin(r) + A*cos(r) + C*sin(r)/r + D*cos(r)/r"""
    rows = []
    rhs_vec = []
    for i in range(len(r_arr)):
        rv = r_arr[i]
        if rv < fs or rv > fe: continue
        u = rv * y_arr[i][ci]
        s, c = sin(rv), cos(rv)
        rows.append([s, c, s/rv, c/rv])
        rhs_vec.append(u)
    n = len(rows)
    MtM = matrix(4, 4)
    Mtb = matrix(4, 1)
    for i in range(n):
        for j in range(4):
            Mtb[j, 0] += rows[i][j] * rhs_vec[i]
            for k in range(4):
                MtM[j, k] += rows[i][j] * rows[i][k]
    x = lu_solve(MtM, Mtb)
    return x[0], x[1], x[2], x[3], n

def extract_6term(r_arr, y_arr, ci, fs, fe):
    """6-term fit: u = B*sin + A*cos + C*sin/r + D*cos/r + E*sin/r² + F*cos/r²"""
    rows = []
    rhs_vec = []
    for i in range(len(r_arr)):
        rv = r_arr[i]
        if rv < fs or rv > fe: continue
        u = rv * y_arr[i][ci]
        s, c = sin(rv), cos(rv)
        rows.append([s, c, s/rv, c/rv, s/rv**2, c/rv**2])
        rhs_vec.append(u)
    n = len(rows)
    MtM = matrix(6, 6)
    Mtb = matrix(6, 1)
    for i in range(n):
        for j in range(6):
            Mtb[j, 0] += rows[i][j] * rhs_vec[i]
            for k in range(6):
                MtM[j, k] += rows[i][j] * rows[i][k]
    x = lu_solve(MtM, Mtb)
    return x[0], x[1], n

print(f"\n{'='*78}")
print("  TAIL AMPLITUDE EXTRACTION")
print("="*78)

ln3 = log(3)
B2_known = mpf(1)/2 - ln3/8
A2_known = pi/8
B3_known = mpf('0.012615939290114712')
A3_known = mpf('-0.215712018451450166')

names = ['f2', 'f3', 'f4', 'f5']
comp_indices = [0, 2, 4, 6]

fit_windows = [
    (500, 2900),
    (800, 2900),
    (1000, 2800),
    (1200, 2600),
]

for ci, name in zip(comp_indices, names):
    print(f"\n  --- {name} ---")
    for fs, fe in fit_windows:
        B4, A4, C4, D4, n4 = extract_4term(r_save, y_save, ci, fs, fe)
        B6, A6, n6 = extract_6term(r_save, y_save, ci, fs, fe)
        print(f"    [{fs},{fe}] 4t: B={nstr(B4,16):>20s}  A={nstr(A4,16):>20s}  (n={n4})")
        print(f"    [{fs},{fe}] 6t: B={nstr(B6,16):>20s}  A={nstr(A6,16):>20s}")

# Best: 6-term fit on wide window
print(f"\n{'='*78}")
print("  BEST VALUES")
print("="*78)

best_fs, best_fe = 800, 2900

amplitudes = {}
for ci, name in zip(comp_indices, names):
    B6, A6, n6 = extract_6term(r_save, y_save, ci, best_fs, best_fe)
    B4, A4, _, _, n4 = extract_4term(r_save, y_save, ci, best_fs, best_fe)
    amplitudes[name] = {'B6': B6, 'A6': A6, 'B4': B4, 'A4': A4}
    print(f"  {name}: B(6t) = {nstr(B6, 18):>22s}   B(4t) = {nstr(B4, 18):>22s}")
    print(f"  {name}: A(6t) = {nstr(A6, 18):>22s}   A(4t) = {nstr(A4, 18):>22s}")

B2_6t = amplitudes['f2']['B6']
A2_6t = amplitudes['f2']['A6']
B3_6t = amplitudes['f3']['B6']
A3_6t = amplitudes['f3']['A6']

print(f"\n  Verification against known values:")
print(f"    B2(6t) - B2_known = {nstr(B2_6t - B2_known, 5)}")
print(f"    A2(6t) - A2_known = {nstr(A2_6t - A2_known, 5)}")
print(f"    B3(6t) - B3_known = {nstr(B3_6t - B3_known, 5)}")
print(f"    A3(6t) - A3_known = {nstr(A3_6t - A3_known, 5)}")
print(f"    B2(4t) - B2_known = {nstr(amplitudes['f2']['B4'] - B2_known, 5)}")
print(f"    A2(4t) - A2_known = {nstr(amplitudes['f2']['A4'] - A2_known, 5)}")
print(f"    B3(4t) - B3_known = {nstr(amplitudes['f3']['B4'] - B3_known, 5)}")
print(f"    A3(4t) - A3_known = {nstr(amplitudes['f3']['A4'] - A3_known, 5)}")

# Use 6-term fit for B4, A4, B5, A5
B4_val = amplitudes['f4']['B6']
A4_val = amplitudes['f4']['A6']
B5_val = amplitudes['f5']['B6']
A5_val = amplitudes['f5']['A6']

# Also compute with 4-term for comparison
B4_4t = amplitudes['f4']['B4']
A4_4t = amplitudes['f4']['A4']
B5_4t = amplitudes['f5']['B4']
A5_4t = amplitudes['f5']['A4']

print(f"\n  B4(6t) = {nstr(B4_val, 18)}")
print(f"  B4(4t) = {nstr(B4_4t, 18)}")
print(f"  A4(6t) = {nstr(A4_val, 18)}")
print(f"  A4(4t) = {nstr(A4_4t, 18)}")
print(f"  B5(6t) = {nstr(B5_val, 18)}")
print(f"  B5(4t) = {nstr(B5_4t, 18)}")
print(f"  A5(6t) = {nstr(A5_val, 18)}")
print(f"  A5(4t) = {nstr(A5_4t, 18)}")

# --- alpha_4, alpha_5 ---
print(f"\n{'='*78}")
print("  ALPHA_4, ALPHA_5")
print("="*78)

A2k = A2_known
B2k = B2_known

# Use both 6t and 4t values
for label, B4v, A4v, B5v, A5v in [
    ("6-term", B4_val, A4_val, B5_val, A5_val),
    ("4-term", B4_4t, A4_4t, B5_4t, A5_4t),
]:
    alpha4 = B4v + A2k*A3_known - B2k*A2k**2/2
    alpha5 = B5v + A3_known**2/2 + A2k*A4v - B2k*A2k*A3_known + A2k**2*(B2k**2 - B3_known)/2 - A2k**4/8
    print(f"\n  [{label}]")
    print(f"    alpha_4 = {nstr(alpha4, 18)}")
    print(f"    alpha_5 = {nstr(alpha5, 18)}")

# --- PSLQ on best values ---
print(f"\n{'='*78}")
print("  PSLQ SEARCHES (using 6-term values)")
print("="*78)

def try_pslq(value, name, max_coeff=1000):
    print(f"\n  --- {name} = {nstr(value, 18)} ---")

    # Rational
    rel = pslq([value, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    Rational: {rel[0]}*x + {rel[1]} = 0  =>  x = {float(mpf(-rel[1])/rel[0]):.15f}")

    # pi
    rel = pslq([value, pi, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    pi: {rel[0]}*x + {rel[1]}*pi + {rel[2]} = 0")

    # pi^2
    rel = pslq([value, pi**2, pi, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    pi^2: {rel[0]}*x + {rel[1]}*pi^2 + {rel[2]}*pi + {rel[3]} = 0")

    # ln2, ln3
    rel = pslq([value, log(2), log(3), mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    ln: {rel[0]}*x + {rel[1]}*ln2 + {rel[2]}*ln3 + {rel[3]} = 0")

    # pi, ln2, ln3
    rel = pslq([value, pi, log(2), log(3), mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    pi+ln: {rel[0]}*x + {rel[1]}*pi + {rel[2]}*ln2 + {rel[3]}*ln3 + {rel[4]} = 0")

    # zeta(3), pi^2, pi, 1
    rel = pslq([value, zeta(3), pi**2, pi, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    zeta3: {rel[0]}*x + {rel[1]}*zeta(3) + {rel[2]}*pi^2 + {rel[3]}*pi + {rel[4]} = 0")

    # Full: pi, pi^2, pi^3, ln2, ln3, zeta3, G, 1
    rel = pslq([value, pi, pi**2, pi**3, log(2), log(3), zeta(3), catalan, mpf(1)], maxcoeff=500)
    if rel:
        labels = ['x', 'pi', 'pi^2', 'pi^3', 'ln2', 'ln3', 'zeta3', 'G', '1']
        terms = ' + '.join(f'{rel[i]}*{labels[i]}' for i in range(len(rel)) if rel[i] != 0)
        print(f"    Full: {terms} = 0")

    # pi, pi^2, ln2, ln3, ln2*ln3, zeta3, 1
    rel = pslq([value, pi, pi**2, log(2), log(3), log(2)*log(3), zeta(3), mpf(1)], maxcoeff=500)
    if rel:
        labels = ['x', 'pi', 'pi^2', 'ln2', 'ln3', 'ln2*ln3', 'zeta3', '1']
        terms = ' + '.join(f'{rel[i]}*{labels[i]}' for i in range(len(rel)) if rel[i] != 0)
        print(f"    Cross: {terms} = 0")

# Best alpha values (average of 6t and 4t as cross-check)
alpha4_6t = B4_val + A2k*A3_known - B2k*A2k**2/2
alpha5_6t = B5_val + A3_known**2/2 + A2k*A4_val - B2k*A2k*A3_known + A2k**2*(B2k**2 - B3_known)/2 - A2k**4/8

try_pslq(B4_val, "B4")
try_pslq(A4_val, "A4")
try_pslq(B5_val, "B5")
try_pslq(A5_val, "A5")
try_pslq(alpha4_6t, "alpha_4")
try_pslq(alpha5_6t, "alpha_5")

# Differences from known f1-only parts
try_pslq(A4_val - 5*pi/128, "A4 - 5pi/128")
try_pslq(A5_val + pi/40, "A5 + pi/40")

# Try B4 - B4_T1 (known from r6_c35: B4_T1 = 0.000728520702...)
B4_T1 = mpf('0.000728520702')
try_pslq(B4_val - B4_T1, "B4 - B4_T1 (f2/f3-dependent part)")

print(f"\n{'='*78}")
print("  FINAL SUMMARY")
print("="*78)
print(f"  B4 = {nstr(B4_val, 16)} (6-term)")
print(f"  B4 = {nstr(B4_4t, 16)} (4-term)")
print(f"  A4 = {nstr(A4_val, 16)} (6-term)")
print(f"  A4 = {nstr(A4_4t, 16)} (4-term)")
print(f"  B5 = {nstr(B5_val, 16)} (6-term)")
print(f"  B5 = {nstr(B5_4t, 16)} (4-term)")
print(f"  A5 = {nstr(A5_val, 16)} (6-term)")
print(f"  A5 = {nstr(A5_4t, 16)} (4-term)")
print(f"  alpha_4 = {nstr(alpha4_6t, 16)}")
print(f"  alpha_5 = {nstr(alpha5_6t, 16)}")
print(f"\n{'='*78}")
