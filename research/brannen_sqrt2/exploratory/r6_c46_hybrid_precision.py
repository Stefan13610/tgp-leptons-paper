#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c46_hybrid_precision.py:
Hybrid approach: float64 ODE (fast, R_max=10000) + mpmath tail fitting.

Key insight: float64 ODE with DOP853, rtol=2.3e-14 gives ~14 digit solutions.
The bottleneck was never ODE precision - it was the tail fitting arithmetic.
By doing the fit at mpmath dps=30, we eliminate the float64 fitting limitation.

Strategy:
1. scipy DOP853 at float64, R_max=10000 (very fast, ~10 seconds)
2. Import tail data into mpmath
3. Multi-term fit at dps=30: 4, 6, 8 terms
4. The wide fitting window [3000, 9000] gives O(1/r⁴) residual ~ 10⁻¹⁴

Expected precision: ~10-12 digits for B4, A4, B5, A5.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math
import time

# --- Phase 1: float64 ODE ---
print("="*78)
print("  PHASE 1: Float64 ODE with scipy DOP853")
print("="*78)

def f1(r):
    if r < 1e-10: return 1.0 - r**2/6 + r**4/120
    return np.sin(r)/r

def f1p(r):
    if r < 1e-10: return -r/3 + r**3/30 - r**5/840
    return np.cos(r)/r - np.sin(r)/r**2

def rhs64(r, y):
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

R_max = 10000
r0 = 1e-6
r_eval = np.arange(1.0, R_max + 0.5, 1.0)  # save at integer r

t0 = time.time()
sol = solve_ivp(rhs64, [r0, R_max], [0.0]*8, method='DOP853',
                rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)
t_ode = time.time() - t0
print(f"\n  {sol.message}, R_max={R_max}")
print(f"  {len(r_eval)} evaluation points, {sol.nfev} function evaluations")
print(f"  Time: {t_ode:.1f}s")

r_data = sol.t
y_data = sol.y  # shape (8, N)

# --- Phase 2: mpmath fitting ---
print(f"\n{'='*78}")
print("  PHASE 2: mpmath multi-term tail fitting (dps=30)")
print("="*78)

from mpmath import mp, mpf, sin as mpsin, cos as mpcos, pi as mppi
from mpmath import log as mplog, matrix, lu_solve, nstr, pslq, zeta, catalan

mp.dps = 30

# Convert data to mpmath
print(f"\n  Converting {len(r_data)} points to mpmath...")
r_mp = [mpf(float(r_data[i])) for i in range(len(r_data))]
# y_mp[i] = [f2, f2', f3, f3', f4, f4', f5, f5'] at r_data[i]
y_mp = []
for i in range(len(r_data)):
    y_mp.append([mpf(float(y_data[j, i])) for j in range(8)])

# --- Multi-term fit with QR via Householder ---
def multiterm_fit(r_arr, y_arr, ci, fs, fe, n_terms):
    """Multi-term fit using normal equations at dps=30."""
    rows = []
    rhs_vec = []
    for i in range(len(r_arr)):
        rv = r_arr[i]
        if rv < fs or rv > fe: continue
        u = rv * y_arr[i][ci]
        s = mpsin(rv)
        c = mpcos(rv)
        row = []
        for k in range(n_terms):
            rk = rv**k
            row.append(s / rk)
            row.append(c / rk)
        rows.append(row)
        rhs_vec.append(u)

    n_data = len(rows)
    n_params = 2 * n_terms

    # Build and solve normal equations
    MtM = matrix(n_params, n_params)
    Mtb = matrix(n_params, 1)
    for i in range(n_data):
        for j in range(n_params):
            Mtb[j, 0] += rows[i][j] * rhs_vec[i]
            for k_col in range(n_params):
                MtM[j, k_col] += rows[i][j] * rows[i][k_col]

    try:
        x = lu_solve(MtM, Mtb)
        return x[0], x[1], n_data
    except:
        return None, None, n_data

# Known values
ln3 = mplog(3)
B2_known = mpf(1)/2 - ln3/8
A2_known = mppi/8
B3_known = mpf('0.012615939290114712')
A3_known = mpf('-0.215712018451450166')

# --- Systematic scan of fit windows and term counts ---
print(f"\n  --- Fit convergence for f2 (B2, A2) ---")
for fs, fe in [(2000, 9000), (3000, 9000), (4000, 9000), (5000, 9000)]:
    print(f"\n    Window [{fs},{fe}]:")
    for nt in [2, 3, 4, 5]:
        B, A, nd = multiterm_fit(r_mp, y_mp, 0, fs, fe, nt)
        if B is None:
            print(f"      {2*nt}t: SINGULAR")
            continue
        be = B - B2_known
        ae = A - A2_known
        print(f"      {2*nt}t: B_err={nstr(be,5):>12s}  A_err={nstr(ae,5):>12s}  (n={nd})")

# --- Main results: scan all components ---
print(f"\n{'='*78}")
print("  RESULTS: all components")
print("="*78)

names = ['f2', 'f3', 'f4', 'f5']
comp_indices = [0, 2, 4, 6]
known_B = [B2_known, B3_known, None, None]
known_A = [A2_known, A3_known, None, None]

best_window = (3000, 9000)
fs, fe = best_window

for ci, name, kb, ka in zip(comp_indices, names, known_B, known_A):
    print(f"\n  === {name} [{fs},{fe}] ===")
    for nt in [2, 3, 4, 5]:
        B, A, nd = multiterm_fit(r_mp, y_mp, ci, fs, fe, nt)
        if B is None:
            print(f"    {2*nt}t: SINGULAR")
            continue
        line = f"    {2*nt}t: B={nstr(B,16):>20s}  A={nstr(A,16):>20s}"
        if kb is not None:
            line += f"  B_err={nstr(B-kb, 5):>12s}"
        if ka is not None:
            line += f"  A_err={nstr(A-ka, 5):>12s}"
        print(line)

# --- Best estimates (8-term fit on [3000, 9000]) ---
print(f"\n{'='*78}")
print("  BEST ESTIMATES")
print("="*78)

best_nt = 4  # 8-parameter fit

results = {}
for ci, name in zip(comp_indices, names):
    B, A, nd = multiterm_fit(r_mp, y_mp, ci, fs, fe, best_nt)
    results[name] = (B, A)
    print(f"  {name}: B = {nstr(B, 18)}, A = {nstr(A, 18)}")

B2_b, A2_b = results['f2']
B3_b, A3_b = results['f3']
B4_b, A4_b = results['f4']
B5_b, A5_b = results['f5']

print(f"\n  Verification:")
print(f"    B2 err = {nstr(B2_b - B2_known, 5)}")
print(f"    A2 err = {nstr(A2_b - A2_known, 5)}")
print(f"    B3 err = {nstr(B3_b - B3_known, 5)}")
print(f"    A3 err = {nstr(A3_b - A3_known, 5)}")

# alpha_4, alpha_5
A2k = A2_known
B2k = B2_known

alpha4 = B4_b + A2k*A3_known - B2k*A2k**2/2
alpha5 = B5_b + A3_known**2/2 + A2k*A4_b - B2k*A2k*A3_known + A2k**2*(B2k**2-B3_known)/2 - A2k**4/8

print(f"\n  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")

# Window stability check
print(f"\n  --- Window stability ({2*best_nt}-term) ---")
for wfs, wfe in [(2000, 9000), (3000, 9000), (4000, 9000), (5000, 9000),
                  (3000, 8000), (3000, 7000), (4000, 8000)]:
    B4w, A4w, _ = multiterm_fit(r_mp, y_mp, 4, wfs, wfe, best_nt)
    B5w, A5w, _ = multiterm_fit(r_mp, y_mp, 6, wfs, wfe, best_nt)
    if B4w is None: continue
    a4w = B4w + A2k*A3_known - B2k*A2k**2/2
    a5w = B5w + A3_known**2/2 + A2k*A4w - B2k*A2k*A3_known + A2k**2*(B2k**2-B3_known)/2 - A2k**4/8
    print(f"    [{wfs},{wfe}]: α4={nstr(a4w,14):>18s}  α5={nstr(a5w,14):>18s}")

# --- PSLQ ---
print(f"\n{'='*78}")
print("  PSLQ SEARCHES")
print("="*78)

def try_pslq(value, name, max_coeff=2000):
    print(f"\n  --- {name} = {nstr(value, 18)} ---")
    for basis_name, basis in [
        ("rational", [value, mpf(1)]),
        ("pi", [value, mppi, mpf(1)]),
        ("pi^2", [value, mppi**2, mppi, mpf(1)]),
        ("pi^2,pi^3", [value, mppi**2, mppi**3, mpf(1)]),
        ("ln2,ln3", [value, mplog(2), mplog(3), mpf(1)]),
        ("pi,ln2,ln3", [value, mppi, mplog(2), mplog(3), mpf(1)]),
        ("zeta3,pi^2,pi", [value, zeta(3), mppi**2, mppi, mpf(1)]),
        ("full8", [value, mppi, mppi**2, mppi**3, mplog(2), mplog(3), zeta(3), catalan, mpf(1)]),
        ("pi*ln", [value, mppi*mplog(2), mppi*mplog(3), mppi, mplog(2), mplog(3), mpf(1)]),
    ]:
        try:
            rel = pslq(basis, maxcoeff=max_coeff)
            if rel:
                print(f"    {basis_name}: {rel}")
        except:
            pass

try_pslq(alpha4, "alpha_4")
try_pslq(alpha5, "alpha_5")
try_pslq(B4_b, "B4")
try_pslq(A4_b, "A4")
try_pslq(B5_b, "B5")
try_pslq(A5_b, "A5")

# Known decompositions
try_pslq(A4_b - 5*mppi/128, "A4 - 5pi/128")
try_pslq(A5_b + mppi/40, "A5 + pi/40")

print(f"\n{'='*78}")
print("  FINAL SUMMARY")
print("="*78)
print(f"  B4 = {nstr(B4_b, 16)}")
print(f"  A4 = {nstr(A4_b, 16)}")
print(f"  B5 = {nstr(B5_b, 16)}")
print(f"  A5 = {nstr(A5_b, 16)}")
print(f"  alpha_4 = {nstr(alpha4, 16)}")
print(f"  alpha_5 = {nstr(alpha5, 16)}")
print(f"\n{'='*78}")
