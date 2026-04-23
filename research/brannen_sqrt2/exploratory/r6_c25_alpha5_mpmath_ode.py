#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c25_alpha5_mpmath_ode.py:
Compute alpha_5 by solving the FULL nonlinear substrate ODE with mpmath.

Strategy:
- Use mpmath at dps=25 for ~20 digit precision
- RK4 with dr=0.001, R_max=300 (~300000 steps)
- Solve for delta = ±{0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20}
- Extract A_tail from tail fit in [100, 280]
- Polynomial fit (eta_sym-1)/d^2 vs d^2 → alpha_3 (verify), alpha_5 (new)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
import time

mp.mp.dps = 25
pi = mp.pi

print("="*78)
print("  ALPHA_5 from full nonlinear ODE (mpmath, dps=25)")
print("="*78)

def rhs(r, g, gp):
    """g'' + (1/g)(g')^2 + (2/r)g' + g = 1"""
    if g < mp.mpf('1e-15'):
        g = mp.mpf('1e-15')
    if r < mp.mpf('1e-15'):
        gpp = (1 - g) / 4
    else:
        gpp = (1 - g) - gp**2/g - 2*gp/r
    return gpp

def solve_ode_mpmath(g0, R_max, dr):
    """Solve ODE from r=eps to R_max with RK4."""
    r = mp.mpf('1e-8')
    g = g0
    gp = mp.mpf(0)

    # Collect tail data for fitting
    r_data = []
    u_data = []  # u = (g-1)*r

    N = int(float(R_max / dr))
    save_interval = max(1, N // 3000)

    for step in range(N):
        # RK4
        k1g = gp
        k1gp = rhs(r, g, gp)

        r2 = r + dr/2
        g2 = g + dr/2*k1g
        gp2 = gp + dr/2*k1gp
        k2g = gp2
        k2gp = rhs(r2, g2, gp2)

        g3 = g + dr/2*k2g
        gp3 = gp + dr/2*k2gp
        k3g = gp3
        k3gp = rhs(r2, g3, gp3)

        r4 = r + dr
        g4 = g + dr*k3g
        gp4 = gp + dr*k3gp
        k4g = gp4
        k4gp = rhs(r4, g4, gp4)

        g = g + dr/6*(k1g + 2*k2g + 2*k3g + k4g)
        gp = gp + dr/6*(k1gp + 2*k2gp + 2*k3gp + k4gp)
        r = r + dr

        # Save data for tail fitting
        if step % save_interval == 0 and float(r) > 80:
            r_data.append(r)
            u_data.append((g - 1) * r)

    return r_data, u_data

def extract_amplitude(r_data, u_data, r_lo, r_hi):
    """Extract A from u ~ A_c*cos(r) + A_s*sin(r) via least squares."""
    # Filter to fitting window
    rr = []
    uu = []
    for ri, ui in zip(r_data, u_data):
        if float(ri) >= r_lo and float(ri) <= r_hi:
            rr.append(ri)
            uu.append(ui)

    n = len(rr)
    if n < 10:
        return None

    # Normal equations for u = c*cos(r) + s*sin(r)
    # [sum(cos^2)  sum(cos*sin)] [c]   [sum(u*cos)]
    # [sum(cos*sin) sum(sin^2) ] [s] = [sum(u*sin)]

    cc = mp.mpf(0); cs = mp.mpf(0); ss = mp.mpf(0)
    uc = mp.mpf(0); us = mp.mpf(0)
    for ri, ui in zip(rr, uu):
        c = mp.cos(ri)
        s = mp.sin(ri)
        cc += c*c
        cs += c*s
        ss += s*s
        uc += ui*c
        us += ui*s

    det = cc*ss - cs*cs
    coef_c = (ss*uc - cs*us) / det
    coef_s = (cc*us - cs*uc) / det

    A = mp.sqrt(coef_c**2 + coef_s**2)
    return A

# Parameters
R_MAX = mp.mpf(300)
dr = mp.mpf('0.002')  # ~150000 steps
deltas = [mp.mpf(d) for d in ['0.05', '0.08', '0.10', '0.12', '0.15', '0.18', '0.20']]

print(f"\n  R_max = {R_MAX}, dr = {dr}, steps = {int(float(R_MAX/dr))}")
print(f"  Deltas: {[str(d) for d in deltas]}")
print(f"  Fitting window: [100, 280]")

alpha3_exact = mp.mpf('0.08972222367362532604749')

results = []
for d in deltas:
    t0 = time.time()

    # Solve +delta and -delta
    rp, up = solve_ode_mpmath(1 + d, R_MAX, dr)
    rm, um = solve_ode_mpmath(1 - d, R_MAX, dr)

    # Extract amplitudes
    Ap = extract_amplitude(rp, up, 100, 280)
    Am = extract_amplitude(rm, um, 100, 280)

    if Ap is None or Am is None:
        print(f"  d={d}: FAILED")
        continue

    eta_p = Ap / d
    eta_m = Am / d
    eta_sym = (eta_p + eta_m) / 2

    ratio = (eta_sym - 1) / d**2
    F = (eta_sym - 1 - alpha3_exact * d**2) / d**4

    dt = time.time() - t0
    print(f"  d={mp.nstr(d,2)}  eta_sym={mp.nstr(eta_sym, 20)}  ratio={mp.nstr(ratio, 18)}  F={mp.nstr(F, 12)}  ({dt:.1f}s)")

    results.append((d, eta_sym, ratio, F))

print(f"\n{'='*78}")
print(f"  POLYNOMIAL FIT")
print(f"{'='*78}")

if len(results) >= 4:
    # Fit (eta_sym-1)/d^2 = a0 + a1*d^2 + a2*d^4 + a3*d^6
    n = len(results)

    # Build matrix (using mpmath for precision)
    for deg in [2, 3]:
        # A * x = b where A[i,j] = d_i^(2j), b[i] = ratio_i
        A = mp.matrix(n, deg+1)
        b = mp.matrix(n, 1)
        for i, (d, eta, ratio, F) in enumerate(results):
            for j in range(deg+1):
                A[i,j] = d**(2*j)
            b[i] = ratio

        # Solve via least squares: (A^T A) x = A^T b
        ATA = A.T * A
        ATb = A.T * b
        coef = mp.lu_solve(ATA, ATb)

        print(f"\n  Degree {deg} fit of (eta_sym-1)/d^2:")
        print(f"    alpha_3 = {mp.nstr(coef[0], 20)}  (err = {mp.nstr(coef[0] - alpha3_exact, 6)})")
        print(f"    alpha_5 = {mp.nstr(coef[1], 20)}")
        if deg >= 2:
            print(f"    alpha_7 = {mp.nstr(coef[2], 15)}")
        if deg >= 3:
            print(f"    alpha_9 = {mp.nstr(coef[3], 12)}")

    # Also fit F(d) = alpha_5 + alpha_7*d^2 + ...
    for deg in [1, 2]:
        A = mp.matrix(n, deg+1)
        b = mp.matrix(n, 1)
        for i, (d, eta, ratio, F) in enumerate(results):
            for j in range(deg+1):
                A[i,j] = d**(2*j)
            b[i] = F

        ATA = A.T * A
        ATb = A.T * b
        coef = mp.lu_solve(ATA, ATb)

        print(f"\n  Degree {deg} fit of F(d):")
        print(f"    alpha_5 = {mp.nstr(coef[0], 20)}")
        if deg >= 1:
            print(f"    alpha_7 = {mp.nstr(coef[1], 15)}")
        if deg >= 2:
            print(f"    alpha_9 = {mp.nstr(coef[2], 12)}")

print("="*78)
