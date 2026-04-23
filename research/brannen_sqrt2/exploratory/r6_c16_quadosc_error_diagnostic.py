#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c16_quadosc_error_diagnostic.py:
Diagnostic: compute A_1 via quadosc at multiple maxdegree settings
to measure quadosc convergence behavior. Also compute A_1 by manual
summation of half-period integrals with Richardson extrapolation,
to independently verify quadosc.
"""
import sys, io, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 40
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)
Ci = mp.ci
Si = mp.si

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360 - s**9/3991680
    return (s*mp.cos(s) - mp.sin(s))/s**2

def Aker(t): return t * fp(t)**2

def Phi1(t):
    return -Ci(t)/2 + Ci(3*t)/2 - (mp.sin(t) + mp.sin(3*t))/(4*t)

def A1_int(t):
    return mp.cos(t) * Aker(t) * Phi1(t)

print("="*72)
print("  QUADOSC ERROR DIAGNOSTIC for A_1 (dps=40)")
print("="*72)

# --- Test 1: Standard quadosc with different maxdegree ---
print("\n  [1] quadosc with various maxdegree:")
for md in [6, 8, 10, 12]:
    t0 = time.time()
    try:
        val = mp.quadosc(A1_int, [0, mp.inf], period=2*pi, maxdegree=md)
        elapsed = time.time() - t0
        print(f"    maxdeg={md:2d}: A_1 = {mp.nstr(val, 30)}  ({elapsed:.1f}s)")
    except Exception as e:
        print(f"    maxdeg={md:2d}: FAILED: {e}")

# --- Test 2: Manual half-period summation ---
print("\n  [2] Manual half-period summation with Richardson:")

def integrate_half_period(k):
    """Integrate A1_int over [k*pi, (k+1)*pi]"""
    a = mp.mpf(k) * pi
    b = mp.mpf(k + 1) * pi
    return mp.quad(A1_int, [a, b])

# Compute partial sums S_N = sum_{k=0}^{N-1} integral on [k*pi, (k+1)*pi]
# Then apply Richardson/Euler-Maclaurin.

t0 = time.time()
terms = []
partial_sum = mp.mpf(0)
N_max = 200

print(f"    Computing {N_max} half-period integrals ...")
for k in range(N_max):
    ik = integrate_half_period(k)
    partial_sum += ik
    if k < 10 or k % 20 == 19:
        print(f"      S_{k+1:>3d} = {mp.nstr(partial_sum, 25)}  (term = {mp.nstr(ik, 10)})")
    terms.append(ik)

elapsed = time.time() - t0
print(f"    (computed in {elapsed:.1f} s)")

# Apply Euler transform (Richardson) on the alternating tail
print(f"\n    S_200 raw = {mp.nstr(partial_sum, 30)}")

# The tail terms (k >= some K) should alternate in sign and decrease.
# Apply Euler-Maclaurin acceleration.
# Use mpmath's nsum on the sequence for comparison:
# Actually, let's use the shanks transform on partial sums.

# Simple Aitken delta-squared for the last 30 partial sums:
S = [sum(terms[:k+1]) for k in range(N_max)]

# Shanks transform / Wynn epsilon algorithm on last 50 partial sums
from functools import reduce
last_n = 80
sn = [mp.mpf(S[N_max - last_n + i]) for i in range(last_n)]

# Wynn epsilon table
def wynn_epsilon(s):
    n = len(s)
    e = [[mp.mpf(0)]*(n+1) for _ in range(n+1)]
    for i in range(n):
        e[i][1] = s[i]
    for k in range(2, n+1):
        for i in range(n - k + 1):
            diff = e[i+1][k-1] - e[i][k-1]
            if abs(diff) < mp.mpf(10)**(-mp.mp.dps):
                e[i][k] = mp.mpf(10)**(mp.mp.dps)  # overflow guard
            else:
                e[i][k] = e[i+1][k-2] + 1/diff
    # Return best estimates (even columns, last row)
    estimates = []
    for k in range(1, n+1, 2):
        estimates.append(e[0][k])
    return estimates

wy = wynn_epsilon(sn)
print(f"\n    Wynn-epsilon accelerated estimates:")
for i, w in enumerate(wy[:8]):
    print(f"      eps_{2*i+1:2d} = {mp.nstr(w, 28)}")

# Compare
print(f"\n    Wynn best    = {mp.nstr(wy[-2] if len(wy)>1 else wy[0], 28)}")
print(f"    quadosc std  = {mp.nstr(mp.quadosc(A1_int, [0, mp.inf], period=2*pi), 28)}")

print("="*72)
