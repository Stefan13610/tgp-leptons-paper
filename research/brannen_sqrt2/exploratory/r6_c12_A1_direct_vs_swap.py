#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c12_A1_direct_vs_swap.py: Verify that A_1 computed via swap (using Phi_1) equals
A_1 via direct nested evaluation. This pins down whether numerical 1.5e-6 discrepancy
in alpha_3 is:
(a) from the integration (swap vs direct should agree to mpmath precision), or
(b) from alpha_3 really not being exactly pi^2/110.

A_1 direct = int_0^inf cos^2(s) f'(s) J_c(s) ds
where J_c(s) = int_0^s cos(t) t (f'(t))^2 dt.

We use: J_c(s) = J_c(inf) - Phi_c(s), where Phi_c(s) = int_s^inf cos(t) t(f')^2 dt.
J_c(inf) = ln(3)/8 - 1/2. And Phi_c is a well-behaved tail function.

A_1 = J_c(inf) · int cos^2(s) f'(s) ds - int cos^2(s) f'(s) Phi_c(s) ds
    = J_c(inf) · Phi_1(0) - int cos^2(s) f'(s) Phi_c(s) ds

where Phi_1(0) = ln(3)/2 - 1 (shown before).
"""
import sys, io
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
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*mp.cos(s) - mp.sin(s))/s**2

def A_ker(t): return t * fp(t)**2

def Phi1(t):
    return -Ci(t)/2 + Ci(3*t)/2 - (mp.sin(t) + mp.sin(3*t))/(4*t)

# Compute J_c(s) via adaptive integration from 0 to s.
# For efficiency, use ODE approach: we'll build a spline first.
# For dps=40, spline errors matter; use direct quad with cache.

import functools

@functools.lru_cache(maxsize=10000)
def Jc_r_cached(r_str):
    r = mp.mpf(r_str)
    if r < mp.mpf('1e-15'): return mp.mpf(0)
    return mp.quad(lambda t: mp.cos(t) * A_ker(t), [0, r])

def Jc_r(r):
    return Jc_r_cached(mp.nstr(r, 30))

# A_1 direct integrand = cos^2(s) * f'(s) * J_c(s)
def A1_direct_int(s):
    return mp.cos(s)**2 * fp(s) * Jc_r(s)

print("="*72)
print(f"  Direct vs swap: A_1 (slow)")
print("="*72)

# Compute A_1 direct over [0, 60] with quadosc (slow due to nested quad)
mp.mp.dps = 25  # reduce precision for speed
print(f"  A_1 direct (slow, nested, dps={mp.mp.dps}):")

def A1_direct_int_fast(s):
    if s < mp.mpf('1e-10'):
        return mp.mpf(0)
    jc = mp.quad(lambda t: mp.cos(t) * A_ker(t), [0, s])
    return mp.cos(s)**2 * fp(s) * jc

# Truncate at L and check convergence
for L in [30, 60, 100]:
    # Use mp.quad with many nodes along [0, L]
    nodes = [mp.mpf(k)*pi/2 for k in range(int(2*L/pi.__float__())+2) if k*pi/2 <= L] + [mp.mpf(L)]
    nodes = sorted(set(nodes))
    val = mp.quad(A1_direct_int_fast, nodes)
    print(f"    L={L:3d}: A_1 direct = {mp.nstr(val, 18)}")

# A_1 via Phi_1 swap (fast, use Phi_1 closed form)
mp.mp.dps = 40
def A1_swap_int(t):
    return mp.cos(t) * A_ker(t) * Phi1(t)

A1_swap = mp.quadosc(A1_swap_int, [0, mp.inf], period=2*pi)
print(f"\n  A_1 via swap (Phi_1 closed, dps=40) = {mp.nstr(A1_swap, 30)}")

# Alternative formulation: A_1 = J_c(inf) * Phi_1(0) - int cos^2(s) f'(s) Phi_c(s) ds
# Phi_c(s) = int_s^inf cos(t) A(t) dt = J_c(inf) - J_c(s)
# Check this: it's equivalent to swap, just rearranged.

Jc_inf = ln3/8 - mp.mpf(1)/2
Phi1_0 = ln3/2 - 1   # verified earlier

print(f"\n  J_c(inf) = ln(3)/8 - 1/2 = {mp.nstr(Jc_inf, 20)}")
print(f"  Phi_1(0) = ln(3)/2 - 1    = {mp.nstr(Phi1_0, 20)}")
print(f"  J_c(inf)*Phi_1(0)         = {mp.nstr(Jc_inf * Phi1_0, 20)}")

# A1_alt: A_1 = J_c(inf)*Phi_1(0) - int cos^2(s) f'(s) Phi_c(s) ds
# where Phi_c(s) = int_s^inf cos(t) A(t) dt.
# This is just cos-variant of Phi_1 with different weighting.
# Actually no, Phi_c here is NOT Phi_1 — Phi_1 is with cos^2(s), Phi_c is just cos(t)A(t).

def Phi_c_direct(s):
    return mp.quadosc(lambda t: mp.cos(t) * A_ker(t), [s, mp.inf], period=2*pi)

# A_1 = J_c(inf)*Phi_1(0) - int cos^2 f' Phi_c ds
# This rearrangement shouldn't give anything new, but verify.

# Alternative useful form:
# Note cos^2(s)f'(s) Phi_c(s): at s -> inf, Phi_c(s) -> 0 but with slow decay.
# cos^2(s)f'(s) ~ cos^3(s)/s - slow oscillatory decay.
# Product ~ cos^3(s)/s * O(sin/s) ~ O(1/s^2), integrable.

print(f"\n{'='*72}")
print(f"  Checking A_1 via J_c(inf)*Phi_1(0) - correction:")
A1_leading = Jc_inf * Phi1_0
print(f"    A_1 leading (J_c*Phi_1) = {mp.nstr(A1_leading, 20)}")
print(f"    A_1 actual (swap)       = {mp.nstr(A1_swap, 20)}")
correction = A1_leading - A1_swap
print(f"    correction difference    = {mp.nstr(correction, 20)}")
print(f"{'='*72}")

# ==== Now verify my formula for alpha_3 via yet another path ====
# alpha_3 = I_sin^2/2 + P_cos = pi^2/128 + (K_c^(I) - 2 K_c^(II))
# With K_c^(II) via ULTRA-PRECISE swap evaluation.
KcI = (ln2 - 1)/6
print(f"\n  K_c^(I)  = (ln2 - 1)/6 = {mp.nstr(KcI, 30)}")

# Full K_c^(II) via all 4 A_i (using Phi_i closed):
def Phi2(t):
    return (Si(3*t) - Si(t))/2 - (mp.cos(t) - mp.cos(3*t))/(4*t)
def Phi3(t):
    return (3*mp.sin(t) - mp.sin(3*t))/(8*t) + 3*(Ci(3*t) - Ci(t))/8 + (mp.cos(3*t) - mp.cos(t))/(8*t**2)
def Phi4(t):
    return (5*mp.cos(t) - mp.cos(3*t))/(8*t) - pi/8 + 5*Si(t)/8 - 3*Si(3*t)/8 - (mp.sin(t) + mp.sin(3*t))/(8*t**2)

A2 = mp.quadosc(lambda t: mp.sin(t) * A_ker(t) * Phi2(t), [0, mp.inf], period=2*pi)
A3 = mp.quadosc(lambda t: mp.cos(t) * A_ker(t) * Phi3(t), [0, mp.inf], period=2*pi)
A4 = mp.quadosc(lambda t: mp.sin(t) * A_ker(t) * Phi4(t), [0, mp.inf], period=2*pi)

KcII = -A1_swap - A2 + A3 - A4
Pcos = KcI - 2*KcII
alpha3 = pi**2/128 + Pcos

print(f"  K_c^(II) = {mp.nstr(KcII, 30)}")
print(f"  P_cos    = K_c^(I) - 2 K_c^(II) = {mp.nstr(Pcos, 30)}")
print(f"  alpha_3  = pi^2/128 + P_cos      = {mp.nstr(alpha3, 30)}")
print(f"  pi^2/110 target                   = {mp.nstr(pi**2/110, 30)}")
print(f"  alpha_3 - pi^2/110              = {mp.nstr(alpha3 - pi**2/110, 8)}")

# Candidate: maybe alpha_3 involves chi_2 or related?
# Try alpha_3 = pi^2/110 + c·X for small correction X.
diff = alpha3 - pi**2/110
print(f"\n  Correction candidates for (alpha_3 - pi^2/110 = {mp.nstr(diff, 8)}):")
# Check various coefficients
ln2sq = ln2**2
ln3sq = ln3**2
ln2ln3 = ln2*ln3
print(f"    diff / (ln2^2)  = {mp.nstr(diff/ln2sq, 15)}")
print(f"    diff / (ln3^2)  = {mp.nstr(diff/ln3sq, 15)}")
print(f"    diff / (ln2·ln3)= {mp.nstr(diff/ln2ln3, 15)}")
print(f"    diff / pi^2     = {mp.nstr(diff/pi**2, 15)}")
print(f"    diff·128        = {mp.nstr(diff*128, 15)}")
print(f"    diff·7040       = {mp.nstr(diff*7040, 15)}")

print("="*72)
