#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c40_B4A4_mpmath_highprec.py:
High-precision B4, A4, B5, A5 using mpmath arithmetic (dps=25).

Strategy:
1. Solve 8-component perturbation ODE (f2..f5) with mpmath RK4, dps=25
2. Extract tail amplitudes via least-squares fit of r*f_n to B*sin(r)+A*cos(r)
3. Compute alpha_4, alpha_5 to 10+ digits
4. Run PSLQ to search for closed forms

The mpmath RK4 at dps=25 breaks the float64 barrier (~5 digits) that limited
previous computations in r6_c30 through r6_c34.

Estimated runtime: ~5-10 min for R_max=500, dr=0.005 (8-component system)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import time
from mpmath import mp, mpf, sin, cos, sqrt, pi, log, matrix, lu_solve, fabs
from mpmath import nstr

mp.dps = 25

print("="*78)
print("  HIGH-PRECISION B4/A4/B5/A5 via mpmath RK4 (dps=25)")
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

# --- RHS: 8-component system [f2, f2', f3, f3', f4, f4', f5, f5'] ---
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
R_max = 500
dr = mpf('0.005')
N_steps = int(R_max / float(dr))

print(f"\n  R_max = {R_max}, dr = {float(dr)}, N_steps = {N_steps}")
print(f"  dps = {mp.dps}")

r = mpf('1e-6')
y = [mpf(0)]*8

t0 = time.time()

# Store solution at integer r values for tail fitting
r_save = []
y_save = []
save_interval = int(1.0 / float(dr))  # save every 1.0 in r

for step in range(N_steps):
    y = rk4_step(rhs, r, y, dr)
    r += dr

    # Save at approximately integer r values
    if (step + 1) % save_interval == 0:
        r_save.append(r)
        y_save.append(list(y))

    # Progress
    if (step + 1) % (N_steps // 10) == 0:
        elapsed = time.time() - t0
        pct = 100*(step+1)/N_steps
        print(f"    {pct:.0f}% done (r={float(r):.1f}), elapsed={elapsed:.1f}s")

t_total = time.time() - t0
print(f"\n  ODE integration complete in {t_total:.1f}s")
print(f"  Saved {len(r_save)} points")

# --- Extract tail amplitudes via least-squares ---
# u_n(r) = r * f_n(r) ~ B_n*sin(r) + A_n*cos(r) for large r
# Fit on [fit_start, fit_end]

def extract_amplitudes(r_arr, y_arr, comp_idx, fit_start, fit_end):
    """
    Extract B, A from r*f_n ~ B*sin(r) + A*cos(r).
    comp_idx: 0=f2, 2=f3, 4=f4, 6=f5
    """
    # Build matrices for least-squares: M @ [B, A]^T = b
    # Using normal equations: M^T M x = M^T b
    S_ss = mpf(0)  # sum sin^2
    S_cc = mpf(0)  # sum cos^2
    S_sc = mpf(0)  # sum sin*cos
    S_sb = mpf(0)  # sum sin*b
    S_cb = mpf(0)  # sum cos*b

    n_pts = 0
    for i in range(len(r_arr)):
        rv = r_arr[i]
        if rv < fit_start or rv > fit_end:
            continue
        fval = y_arr[i][comp_idx]
        u = rv * fval  # tail variable
        s = sin(rv)
        c = cos(rv)
        S_ss += s*s
        S_cc += c*c
        S_sc += s*c
        S_sb += s*u
        S_cb += c*u
        n_pts += 1

    # Solve 2x2 normal equations
    det = S_ss*S_cc - S_sc**2
    B = (S_cc*S_sb - S_sc*S_cb) / det
    A = (S_ss*S_cb - S_sc*S_sb) / det

    return B, A, n_pts

# Also do 4-term fit: B*sin + A*cos + C*sin/r + D*cos/r
def extract_amplitudes_4term(r_arr, y_arr, comp_idx, fit_start, fit_end):
    """4-term fit: u = B*sin(r) + A*cos(r) + C*sin(r)/r + D*cos(r)/r"""
    rows = []
    rhs_vec = []
    for i in range(len(r_arr)):
        rv = r_arr[i]
        if rv < fit_start or rv > fit_end:
            continue
        fval = y_arr[i][comp_idx]
        u = rv * fval
        s = sin(rv)
        c = cos(rv)
        rows.append([s, c, s/rv, c/rv])
        rhs_vec.append(u)

    n = len(rows)
    # Normal equations: M^T M x = M^T b
    MtM = matrix(4, 4)
    Mtb = matrix(4, 1)
    for i in range(n):
        for j in range(4):
            Mtb[j, 0] += rows[i][j] * rhs_vec[i]
            for k in range(4):
                MtM[j, k] += rows[i][j] * rows[i][k]

    x = lu_solve(MtM, Mtb)
    return x[0], x[1], x[2], x[3], n

print(f"\n{'='*78}")
print("  TAIL AMPLITUDE EXTRACTION")
print("="*78)

# Known values for verification
ln3 = log(3)
B2_known = mpf(1)/2 - ln3/8
A2_known = pi/8

# Multiple fit windows for stability check
fit_windows = [
    (200, 480),
    (250, 480),
    (300, 480),
    (150, 450),
]

names = ['f2', 'f3', 'f4', 'f5']
comp_indices = [0, 2, 4, 6]

for ci, name in zip(comp_indices, names):
    print(f"\n  --- {name} ---")
    for fs, fe in fit_windows:
        B2t, A2t, npts = extract_amplitudes(r_save, y_save, ci, fs, fe)
        B4t, A4t, C4t, D4t, npts4 = extract_amplitudes_4term(r_save, y_save, ci, fs, fe)
        print(f"    [{fs},{fe}] 2-term: B={nstr(B2t,18):>22s}  A={nstr(A2t,18):>22s}  (n={npts})")
        print(f"    [{fs},{fe}] 4-term: B={nstr(B4t,18):>22s}  A={nstr(A4t,18):>22s}  (n={npts4})")

# Best values: use widest 4-term fit
print(f"\n{'='*78}")
print("  BEST VALUES (4-term fit, [200,480])")
print("="*78)

amplitudes = {}
for ci, name in zip(comp_indices, names):
    B, A, C, D, n = extract_amplitudes_4term(r_save, y_save, ci, 200, 480)
    amplitudes[name] = (B, A)
    print(f"  {name}: B = {nstr(B, 20)}")
    print(f"  {name}: A = {nstr(A, 20)}")

B2, A2 = amplitudes['f2']
B3, A3 = amplitudes['f3']
B4, A4 = amplitudes['f4']
B5, A5 = amplitudes['f5']

# Verify against known
print(f"\n  Verification:")
print(f"    B2 = {nstr(B2, 18)}")
print(f"    B2_known = {nstr(B2_known, 18)}")
print(f"    B2 error = {nstr(B2 - B2_known, 5)}")
print(f"    A2 = {nstr(A2, 18)}")
print(f"    A2_known = {nstr(A2_known, 18)}")
print(f"    A2 error = {nstr(A2 - A2_known, 5)}")

# Known B3, A3 (30 and 25 digits respectively)
B3_known = mpf('0.012615939290114712')
A3_known = mpf('-0.215712018451450166')
print(f"    B3 = {nstr(B3, 18)}")
print(f"    B3_known = {nstr(B3_known, 18)}")
print(f"    B3 error = {nstr(B3 - B3_known, 5)}")
print(f"    A3 = {nstr(A3, 18)}")
print(f"    A3_known = {nstr(A3_known, 18)}")
print(f"    A3 error = {nstr(A3 - A3_known, 5)}")

# --- Compute alpha_4, alpha_5 ---
print(f"\n{'='*78}")
print("  ALPHA_4, ALPHA_5")
print("="*78)

# Use known A2, B2 for highest precision
A2k = A2_known
B2k = B2_known

# alpha_4 = B4 + A2*A3 - B2*A2^2/2
alpha4 = B4 + A2k*A3_known - B2k*A2k**2/2
print(f"  alpha_4 = B4 + A2*A3 - B2*A2^2/2")
print(f"  alpha_4 = {nstr(alpha4, 18)}")

# alpha_5 = B5 + A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
alpha5 = B5 + A3_known**2/2 + A2k*A4 - B2k*A2k*A3_known + A2k**2*(B2k**2 - B3_known)/2 - A2k**4/8
print(f"  alpha_5 = B5 + A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8")
print(f"  alpha_5 = {nstr(alpha5, 18)}")

print(f"\n  Previous estimates (float64):")
print(f"    alpha_4 ~ -0.02460")
print(f"    alpha_5 ~ 0.02751")

# --- PSLQ searches ---
print(f"\n{'='*78}")
print("  PSLQ SEARCHES")
print("="*78)

from mpmath import identify, pslq, euler as euler_gamma, catalan, zeta

# Constant basis for PSLQ
constants = {
    'pi': pi,
    'pi^2': pi**2,
    'ln2': log(2),
    'ln3': log(3),
    'ln5': log(5),
    'zeta3': zeta(3),
    'G': catalan,
    'gamma': euler_gamma,
    'pi*ln2': pi*log(2),
    'pi*ln3': pi*log(3),
    'pi^3': pi**3,
    'ln2^2': log(2)**2,
    'ln3^2': log(3)**2,
}

def try_pslq(value, name, max_coeff=1000):
    """Try PSLQ with various constant bases."""
    print(f"\n  --- PSLQ for {name} = {nstr(value, 18)} ---")

    # Try mpmath identify first
    try:
        result = identify(value, tol=1e-12)
        if result:
            print(f"    identify() -> {result}")
    except:
        pass

    # Try rational
    rel = pslq([value, mpf(1)], maxcoeff=max_coeff)
    if rel:
        # rel[0]*value + rel[1] = 0 => value = -rel[1]/rel[0]
        print(f"    Rational: {rel[0]}*x + {rel[1]} = 0 => x = {mpf(-rel[1])/rel[0]}")

    # Try pi relation
    rel = pslq([value, pi, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    pi: {rel[0]}*x + {rel[1]}*pi + {rel[2]} = 0")

    # Try pi^2 relation
    rel = pslq([value, pi**2, pi, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    pi^2: {rel[0]}*x + {rel[1]}*pi^2 + {rel[2]}*pi + {rel[3]} = 0")

    # Try ln2, ln3
    rel = pslq([value, log(2), log(3), mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    ln: {rel[0]}*x + {rel[1]}*ln2 + {rel[2]}*ln3 + {rel[3]} = 0")

    # Try pi, ln2, ln3
    rel = pslq([value, pi, log(2), log(3), mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    pi+ln: {rel[0]}*x + {rel[1]}*pi + {rel[2]}*ln2 + {rel[3]}*ln3 + {rel[4]} = 0")

    # Try zeta(3), pi^2, 1
    rel = pslq([value, zeta(3), pi**2, mpf(1)], maxcoeff=max_coeff)
    if rel:
        print(f"    zeta3: {rel[0]}*x + {rel[1]}*zeta(3) + {rel[2]}*pi^2 + {rel[3]} = 0")

    # Big basis: value, pi, pi^2, ln2, ln3, zeta3, G, 1
    rel = pslq([value, pi, pi**2, log(2), log(3), zeta(3), catalan, mpf(1)], maxcoeff=500)
    if rel:
        print(f"    Full: {rel[0]}*x + {rel[1]}*pi + {rel[2]}*pi^2 + {rel[3]}*ln2 + {rel[4]}*ln3 + {rel[5]}*zeta3 + {rel[6]}*G + {rel[7]} = 0")

# PSLQ on B4, A4, B5, A5, alpha4, alpha5
try_pslq(B4, "B4")
try_pslq(A4, "A4")
try_pslq(B5, "B5")
try_pslq(A5, "A5")
try_pslq(alpha4, "alpha_4")
try_pslq(alpha5, "alpha_5")

# Also try A4 + 5*pi/128 (the f1-only part is known)
A4_rest = A4 - 5*pi/128  # A4 minus the f1-only contribution
try_pslq(A4_rest, "A4_rest (A4 - 5pi/128)")

# A5 + pi/40 (f1-only part)
A5_rest = A5 + pi/40  # A5 minus the f1-only contribution (which is -pi/40)
try_pslq(A5_rest, "A5_rest (A5 + pi/40)")

print(f"\n{'='*78}")
print("  SUMMARY")
print("="*78)
print(f"  B2 = {nstr(B2, 18)}  (known: {nstr(B2_known, 18)})")
print(f"  A2 = {nstr(A2, 18)}  (known: {nstr(A2_known, 18)})")
print(f"  B3 = {nstr(B3, 18)}  (known: {nstr(B3_known, 18)})")
print(f"  A3 = {nstr(A3, 18)}  (known: {nstr(A3_known, 18)})")
print(f"  B4 = {nstr(B4, 18)}")
print(f"  A4 = {nstr(A4, 18)}")
print(f"  B5 = {nstr(B5, 18)}")
print(f"  A5 = {nstr(A5, 18)}")
print(f"  alpha_4 = {nstr(alpha4, 18)}")
print(f"  alpha_5 = {nstr(alpha5, 18)}")
print(f"\n{'='*78}")
