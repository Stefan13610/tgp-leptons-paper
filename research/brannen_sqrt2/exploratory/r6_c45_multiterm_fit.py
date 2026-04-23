#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c45_multiterm_fit.py:
High-precision tail amplitudes via multi-term fitting at mpmath precision.

Key insight: the tail expansion is
  u_n(r) = r*f_n(r) = B*sin(r) + A*cos(r)
           + C₁*sin(r)/r + D₁*cos(r)/r
           + C₂*sin(r)/r² + D₂*cos(r)/r²
           + C₃*sin(r)/r³ + D₃*cos(r)/r³ + ...

With float64, fits beyond 4-6 terms are ill-conditioned.
With mpmath at dps=30, we can fit 8-12 terms stably because
the condition number (~10^9 for 8 terms) is well within
the 30-digit arithmetic precision.

Strategy:
1. mpmath RK4, dps=30, dr=0.02, R_max=3000
2. Multi-term fit at dps=30 with normal equations
3. Check convergence: 4, 6, 8, 10 term fits
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import time
from mpmath import mp, mpf, sin, cos, pi, log, matrix, lu_solve, nstr, pslq, zeta, catalan

mp.dps = 30

print("="*78)
print("  MULTI-TERM TAIL FIT at mpmath dps=30")
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

# --- ODE Integration ---
R_max = 3000
dr = mpf('0.02')
N_steps = int(R_max / float(dr))

print(f"\n  R_max = {R_max}, dr = {float(dr)}, N_steps = {N_steps}")

r_save = []
y_save = []
save_interval = int(1.0 / float(dr))

r = mpf('1e-6')
y = [mpf(0)]*8

t0 = time.time()

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
print(f"  ODE: {t_total:.1f}s, {len(r_save)} points saved")

# --- Multi-term fitting ---
def multiterm_fit(r_arr, y_arr, ci, fs, fe, n_terms):
    """
    Fit u = r*f_n to:
      sum_{k=0}^{n_terms-1} [a_k * sin(r)/r^k + b_k * cos(r)/r^k]

    n_terms basis functions: 2*n_terms parameters
    k=0: B*sin(r) + A*cos(r)  (leading asymptotic)
    k=1: C*sin(r)/r + D*cos(r)/r
    etc.

    Returns B (coeff of sin(r)), A (coeff of cos(r)).
    """
    rows = []
    rhs_vec = []

    for i in range(len(r_arr)):
        rv = r_arr[i]
        if rv < fs or rv > fe: continue
        u = rv * y_arr[i][ci]
        s = sin(rv)
        c = cos(rv)
        row = []
        for k in range(n_terms):
            rk = rv**k
            row.append(s / rk)
            row.append(c / rk)
        rows.append(row)
        rhs_vec.append(u)

    n_data = len(rows)
    n_params = 2 * n_terms

    # Normal equations: (M^T M) x = M^T b
    MtM = matrix(n_params, n_params)
    Mtb = matrix(n_params, 1)

    for i in range(n_data):
        for j in range(n_params):
            Mtb[j, 0] += rows[i][j] * rhs_vec[i]
            for k_col in range(n_params):
                MtM[j, k_col] += rows[i][j] * rows[i][k_col]

    x = lu_solve(MtM, Mtb)
    B = x[0]  # sin(r) coefficient
    A = x[1]  # cos(r) coefficient

    return B, A, n_data

# --- Results ---
print(f"\n{'='*78}")
print("  MULTI-TERM FIT RESULTS")
print("="*78)

ln3 = log(3)
B2_known = mpf(1)/2 - ln3/8
A2_known = pi/8
B3_known = mpf('0.012615939290114712')
A3_known = mpf('-0.215712018451450166')

names = ['f2', 'f3', 'f4', 'f5']
comp_indices = [0, 2, 4, 6]
known_B = [B2_known, B3_known, None, None]
known_A = [A2_known, A3_known, None, None]

fit_window = (800, 2900)
fs, fe = fit_window

for ci, name, kb, ka in zip(comp_indices, names, known_B, known_A):
    print(f"\n  === {name} (fit [{fs},{fe}]) ===")
    for nt in [2, 3, 4, 5, 6, 7, 8]:
        B, A, nd = multiterm_fit(r_save, y_save, ci, fs, fe, nt)
        label = f"{2*nt}t"
        line = f"    {label:>4s}: B={nstr(B,16):>20s}  A={nstr(A,16):>20s}"
        if kb is not None:
            line += f"  B_err={nstr(B-kb, 4):>10s}"
        if ka is not None:
            line += f"  A_err={nstr(A-ka, 4):>10s}"
        print(line)

# --- Best values: highest stable term count ---
# Check which n_terms gives stable B, A by looking at convergence
print(f"\n{'='*78}")
print("  CONVERGENCE ANALYSIS")
print("="*78)

# For f4, f5: track B, A across term counts
for ci, name in [(4, 'f4'), (6, 'f5')]:
    print(f"\n  --- {name} ---")
    prev_B, prev_A = None, None
    for nt in range(2, 9):
        B, A, nd = multiterm_fit(r_save, y_save, ci, fs, fe, nt)
        if prev_B is not None:
            dB = B - prev_B
            dA = A - prev_A
            print(f"    {2*nt:2d}t: B={nstr(B,16):>20s}  ΔB={nstr(dB,4):>10s}  A={nstr(A,16):>20s}  ΔA={nstr(dA,4):>10s}")
        else:
            print(f"    {2*nt:2d}t: B={nstr(B,16):>20s}  {'':>10s}  A={nstr(A,16):>20s}")
        prev_B, prev_A = B, A

# --- Best estimates ---
print(f"\n{'='*78}")
print("  BEST ESTIMATES")
print("="*78)

# Use the highest term count that's stable
best_nt = 6  # 12-parameter fit

B2_best, A2_best, _ = multiterm_fit(r_save, y_save, 0, fs, fe, best_nt)
B3_best, A3_best, _ = multiterm_fit(r_save, y_save, 2, fs, fe, best_nt)
B4_best, A4_best, _ = multiterm_fit(r_save, y_save, 4, fs, fe, best_nt)
B5_best, A5_best, _ = multiterm_fit(r_save, y_save, 6, fs, fe, best_nt)

print(f"  Using {2*best_nt}-term fit on [{fs},{fe}]:")
print(f"  B2 = {nstr(B2_best, 18)}  (err: {nstr(B2_best - B2_known, 5)})")
print(f"  A2 = {nstr(A2_best, 18)}  (err: {nstr(A2_best - A2_known, 5)})")
print(f"  B3 = {nstr(B3_best, 18)}  (err: {nstr(B3_best - B3_known, 5)})")
print(f"  A3 = {nstr(A3_best, 18)}  (err: {nstr(A3_best - A3_known, 5)})")
print(f"  B4 = {nstr(B4_best, 18)}")
print(f"  A4 = {nstr(A4_best, 18)}")
print(f"  B5 = {nstr(B5_best, 18)}")
print(f"  A5 = {nstr(A5_best, 18)}")

# alpha_4, alpha_5
A2k = A2_known
B2k = B2_known

alpha4 = B4_best + A2k*A3_known - B2k*A2k**2/2
alpha5 = B5_best + A3_known**2/2 + A2k*A4_best - B2k*A2k*A3_known + A2k**2*(B2k**2-B3_known)/2 - A2k**4/8

print(f"\n  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")

# Also try multiple fit windows for stability
print(f"\n  --- Fit window stability (using {2*best_nt}-term) ---")
windows = [(500, 2900), (800, 2900), (1000, 2800), (600, 2500), (800, 2500)]
for wfs, wfe in windows:
    B4w, A4w, _ = multiterm_fit(r_save, y_save, 4, wfs, wfe, best_nt)
    B5w, A5w, _ = multiterm_fit(r_save, y_save, 6, wfs, wfe, best_nt)
    a4w = B4w + A2k*A3_known - B2k*A2k**2/2
    a5w = B5w + A3_known**2/2 + A2k*A4w - B2k*A2k*A3_known + A2k**2*(B2k**2-B3_known)/2 - A2k**4/8
    print(f"    [{wfs},{wfe}]: B4={nstr(B4w,14):>18s} A4={nstr(A4w,14):>18s} α4={nstr(a4w,12):>16s} α5={nstr(a5w,12):>16s}")

# --- PSLQ on best values ---
print(f"\n{'='*78}")
print("  PSLQ")
print("="*78)

def try_pslq(value, name, max_coeff=1000):
    print(f"\n  --- {name} = {nstr(value, 18)} ---")
    for basis_name, basis in [
        ("pi", [value, pi, mpf(1)]),
        ("pi^2", [value, pi**2, pi, mpf(1)]),
        ("ln2,ln3", [value, log(2), log(3), mpf(1)]),
        ("pi,ln2,ln3", [value, pi, log(2), log(3), mpf(1)]),
        ("zeta3,pi^2,pi", [value, zeta(3), pi**2, pi, mpf(1)]),
        ("full", [value, pi, pi**2, pi**3, log(2), log(3), zeta(3), catalan, mpf(1)]),
    ]:
        try:
            rel = pslq(basis, maxcoeff=max_coeff)
            if rel:
                print(f"    {basis_name}: {rel}")
        except:
            pass

try_pslq(alpha4, "alpha_4")
try_pslq(alpha5, "alpha_5")
try_pslq(B4_best, "B4")
try_pslq(A4_best, "A4")
try_pslq(A4_best - 5*pi/128, "A4 - 5pi/128")

print(f"\n{'='*78}")
