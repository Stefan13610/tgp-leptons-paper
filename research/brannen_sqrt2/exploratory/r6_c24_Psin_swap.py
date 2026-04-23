#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c24_Psin_swap.py:
Compute P_sin = K_s^(I) - 2*K_s^(II) via the Fubini swap.

Recap:
  P_cos = K_c^(I) - 2*K_c^(II) = (ln2-1)/6 - 2*K_c^(II)
  P_sin = K_s^(I) - 2*K_s^(II) = pi/12 - 2*K_s^(II)

where K_s^(I) = pi/12 (PROVEN: confirmed to 32 digits)

K_s^(II) = int_0^inf sin(s) * s * f1'(s) * h'(s) ds

Using the same Fubini swap as for K_c^(II) (r6_c6), but with sin(s):
h'(s) = -u2'(s)/s + u2(s)/s^2  (where u2 = -r*h, sign convention may differ)

Actually, let me re-derive from the swap formula:
u2(r) = int_0^r sin(r-t) * source2(t) dt  where source2(t) = -t*(f1'(t))^2
u2'(r) = int_0^r cos(r-t) * source2(t) dt

f2 = u2/r, f2' = (u2' - u2/r)/r = (r*u2' - u2)/r^2

K_s^(II) = int_0^inf sin(s)*s*f1'(s)*f2'(s) ds
         = int_0^inf sin(s)*f1'(s)*(s*u2'(s) - u2(s))/s ds
         = int_0^inf sin(s)*f1'(s)*[u2'(s) - u2(s)/s] ds

Substituting u2(s) and u2'(s):
u2(s) = int_0^s sin(s-t)*src2(t) dt = sin(s)*C2(s) - cos(s)*S2(s)
u2'(s) = cos(s)*C2(s) + sin(s)*S2(s)

u2'(s) - u2(s)/s = cos(s)*C2(s) + sin(s)*S2(s) - [sin(s)*C2(s) - cos(s)*S2(s)]/s
                  = C2(s)*[cos(s) - sin(s)/s] + S2(s)*[sin(s) + cos(s)/s]

where C2(s) = int_0^s cos(t)*src2(t) dt, S2(s) = int_0^s sin(t)*src2(t) dt

So K_s^(II) = int_0^inf sin(s)*f1'(s) * {C2(s)*[cos(s)-sin(s)/s] + S2(s)*[sin(s)+cos(s)/s]} ds

Swapping (Fubini), using src2(t) = -t*(f1'(t))^2:
K_s^(II) = int_0^inf src2(t) * [cos(t)*Psi_1(t) + sin(t)*Psi_2(t)] dt

Wait, let me be more careful. Let me compute this via the same decomposition
as K_c^(II), which was done in r6_c6_KcII_phi_swap.py.

The key decomposition for K_c^(II) was:
K_c^(II) = -A1 - A2 + A3_kc - A4
where A_i = int_0^inf tau_i(t) * t*(f1'(t))^2 * Phi_i(t) dt

For K_s^(II), I need the SAME structure but with sin(s) instead of cos(s)
in the outer integral. This changes the tail functions from Phi_i to Psi_i.

Actually, let me just compute K_s^(II) directly via quadosc.
I need f2'(s) at each s, which requires evaluating the nested integral C2(s), S2(s).

Alternative: compute via the A_i swap decomposition.

The MOST DIRECT approach: use mpmath quadosc with f2'(s) computed via quadosc
at each evaluation point. This is a nested quadrature.

Let me try it.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 35

pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

def f1(r):
    if r < mp.mpf('1e-30'): return mp.mpf(1)
    return mp.sin(r) / r

def f1p(r):
    if r < mp.mpf('1e-30'): return mp.mpf(0)
    return mp.cos(r)/r - mp.sin(r)/r**2

def src2(t):
    """Source for u2: -t*(f1')^2"""
    return -t * f1p(t)**2

# C2(s) = int_0^s cos(t)*src2(t) dt
# S2(s) = int_0^s sin(t)*src2(t) dt
# u2(s) = sin(s)*C2(s) - cos(s)*S2(s)
# u2'(s) = cos(s)*C2(s) + sin(s)*S2(s)
# f2(s) = u2(s)/s
# f2'(s) = (u2'(s) - u2(s)/s)/1 = [r*u2'(s) - u2(s)]/r^2 ... hmm wait
# f2'(s) = (u2'(s)*s - u2(s))/s^2

# For the swap decomposition of K_s^(II):
# K_s^(II) = int_0^inf sin(s)*s*f1'(s)*f2'(s) ds
#          = int_0^inf sin(s)*f1'(s) * (s*f2'(s)) ds
#          = int_0^inf sin(s)*f1'(s) * (u2'(s) - u2(s)/s) ds

# Substituting u2, u2':
# = int_0^inf sin(s)*f1'(s) * [cos(s)*C2(s) + sin(s)*S2(s)
#                              - (sin(s)*C2(s) - cos(s)*S2(s))/s] ds

# = int_0^inf sin(s)*f1'(s) * [C2(s)*(cos(s) - sin(s)/s) + S2(s)*(sin(s) + cos(s)/s)] ds

# Using C2(s) = int_0^s cos(t)*src2(t)dt, swap order:
# = int_0^inf src2(t) * cos(t) * [int_t^inf sin(s)*f1'(s)*(cos(s)-sin(s)/s) ds] dt
# + int_0^inf src2(t) * sin(t) * [int_t^inf sin(s)*f1'(s)*(sin(s)+cos(s)/s) ds] dt

# With src2(t) = -t*(f1'(t))^2, this becomes:
# K_s^(II) = -int_0^inf t*(f1'(t))^2 * [cos(t)*Psi_1(t) + sin(t)*Psi_2(t)] dt
# where:
# Psi_1(t) = int_t^inf sin(s)*f1'(s)*(cos(s) - sin(s)/s) ds
# Psi_2(t) = int_t^inf sin(s)*f1'(s)*(sin(s) + cos(s)/s) ds

# Let me compute the Psi functions analytically (or at least reduce them).

# sin(s)*f1'(s) = sin(s)*[cos(s)/s - sin(s)/s^2]
#               = sin(s)*cos(s)/s - sin^2(s)/s^2

# Psi_1(t):
# sin(s)*f1'(s)*(cos(s) - sin(s)/s)
# = [sin(s)cos(s)/s - sin^2(s)/s^2] * [cos(s) - sin(s)/s]
# = sin(s)cos^2(s)/s - sin^2(s)cos(s)/s^2 - sin(s)cos(s)sin(s)/s^2 + sin^3(s)/s^3
# = sin(s)cos^2(s)/s - sin^2(s)cos(s)/s^2 - sin^2(s)cos(s)/s^2 + sin^3(s)/s^3

# Hmm, this is getting messy. Let me instead take a hybrid approach:
# 1. Compute K_s^(II) numerically using direct nested quadrature (slow but reliable)
# 2. Verify with a swap approach

print("="*78)
print("  Computing P_sin via direct nested quadrature")
print("="*78)

# Approach 1: Direct computation using mpmath ODE for u2 and then quadosc

# First, build u2(s) and u2'(s) via their Green's function integrals.
# Cache the partial integrals.

from functools import lru_cache

# For evaluating u2(s) at each s, I use the Green's function:
# u2(s) = sin(s)*C2(s) - cos(s)*S2(s)
# I need C2(s) = int_0^s cos(t)*src2(t) dt

# Direct quadrature for C2(s) at each s is expensive but doable.

def u2_and_u2p(s):
    """Compute u2(s) and u2'(s) using Green's function integrals."""
    if s < mp.mpf('1e-10'):
        return mp.mpf(0), mp.mpf(0)
    C2 = mp.quad(lambda t: mp.cos(t)*src2(t), [0, s])
    S2 = mp.quad(lambda t: mp.sin(t)*src2(t), [0, s])
    u2 = mp.sin(s)*C2 - mp.cos(s)*S2
    u2p = mp.cos(s)*C2 + mp.sin(s)*S2
    return u2, u2p

def f2p(s):
    """f2'(s) = (u2'*s - u2)/s^2"""
    if s < mp.mpf('1e-10'):
        return mp.mpf(0)
    u2, u2p = u2_and_u2p(s)
    return (u2p * s - u2) / s**2

# K_s^(II) = int_0^inf sin(s)*s*f1'(s)*f2'(s) ds
# This is expensive because f2'(s) requires a quadrature at each s.
# Let me test with a finite upper limit first.

print("\n  Step 1: Verify K_c^(II) via direct nested quadrature...")
print("  (This confirms our method works correctly)")

# K_c^(II) = int_0^inf cos(s)*s*f1'(s)*f2'(s) ds
def KcII_integrand(s):
    return mp.cos(s) * s * f1p(s) * f2p(s)

KcII_known = mp.mpf('-0.031879037931728580134156')

# Test at finite limits to see convergence
for Rmax in [20, 50]:
    val = mp.quad(KcII_integrand, [0, Rmax], maxdegree=7)
    print(f"    int_0^{Rmax}: {mp.nstr(val, 20)} (diff from known: {mp.nstr(val - KcII_known, 6)})")

# Now compute K_s^(II)
print("\n  Step 2: Compute K_s^(II) = int sin(s)*s*f1'(s)*f2'(s) ds ...")

def KsII_integrand(s):
    return mp.sin(s) * s * f1p(s) * f2p(s)

for Rmax in [20, 50]:
    val = mp.quad(KsII_integrand, [0, Rmax], maxdegree=7)
    print(f"    int_0^{Rmax}: {mp.nstr(val, 20)}")

# Try quadosc for the full integral
print("\n  Step 3: quadosc for K_s^(II)...")
try:
    KsII_full = mp.quadosc(KsII_integrand, [0, mp.inf], omega=1)
    print(f"    K_s^(II) = {mp.nstr(KsII_full, 25)}")
except Exception as e:
    print(f"    quadosc failed: {e}")
    # Fall back to finite limit with extrapolation
    vals = []
    for Rmax in [30, 50, 80, 120]:
        val = mp.quad(KsII_integrand, [0, Rmax], maxdegree=8)
        vals.append((Rmax, val))
        print(f"    int_0^{Rmax}: {mp.nstr(val, 20)}")
    KsII_full = vals[-1][1]

# K_s^(I) = pi/12 (proven)
KsI = pi / 12

# P_sin = K_s^(I) - 2*K_s^(II)
Psin = KsI - 2*KsII_full
print(f"\n  K_s^(I) = pi/12 = {mp.nstr(KsI, 25)}")
print(f"  K_s^(II) = {mp.nstr(KsII_full, 25)}")
print(f"  P_sin = K_s^(I) - 2*K_s^(II) = {mp.nstr(Psin, 25)}")

# From the convention: A3 = -P_sin (cos amplitude of u3's tail, negated for sin integral)
# Wait, actually: u3 ~ sin(r)*B3 + cos(r)*(-B3_sin)
# where B3 = P_cos (sin amplitude) and A3 = -Js3 where Js3 = int sin(s)*src3(s) ds
# = P_sin... hmm, I need to be careful.

# From variation of parameters:
# u3(r) ~ sin(r) * int_0^inf cos(s)*src3(s) ds - cos(r) * int_0^inf sin(s)*src3(s) ds
# The sin-amplitude of u3 = int cos(s)*src3 ds = -2*K_c^(II) + K_c^(I) = P_cos = B3
# The cos-amplitude of u3 = -int sin(s)*src3 ds = -(−2*K_s^(II) + K_s^(I)) = 2*K_s^(II) - K_s^(I) = -P_sin
# So: u3 ~ B3*sin(r) + (-P_sin)*cos(r) = B3*sin(r) - P_sin*cos(r)
# In our notation: A3 = -P_sin (cos-amplitude of u3)

# Wait, I said A3 = -Js3 where Js3 = int sin(s)*src3(s) ds = P_sin
# but the cos coefficient is -Js3 = -P_sin
# So A3 (cos-amp of u3) = -P_sin

# Hmm, that gives A3 = -P_sin and in the alpha_5 formula I used A3.
# Let me verify: from r6_c22, the formula for alpha_5 uses:
# p = A2, q = A3, r = A4
# where A_n is the cos-amplitude of u_n.

# A2 = -Js2 = -I_sin = -(-pi/8) = pi/8 ✓
# A3 = -Js3 = -(P_sin) = -P_sin
# But P_sin = K_s^(I) - 2*K_s^(II)

A3_val = -Psin
print(f"\n  A3 = -P_sin = {mp.nstr(A3_val, 25)}")

# Now I need B4, A4 from the f4 equation.
# src4 = -r*[(f2')^2 + 2*f1'*f3' - 2*f1*f1'*f2' - f2*(f1')^2 + f1^2*(f1')^2]

# B4 = int cos(s)*src4(s) ds (sin-amplitude of u4)
# A4 = -int sin(s)*src4(s) ds (cos-amplitude of u4)

# The source has 5 terms:
# (1) (f2')^2      → requires f2'  (single swap level)
# (2) 2*f1'*f3'    → requires f3'  (double swap level - HARD)
# (3) -2*f1*f1'*f2' → requires f2'  (single swap level)
# (4) -f2*(f1')^2  → requires f2   (single swap level)
# (5) f1^2*(f1')^2 → pure f1 (elementary)

# Let me compute each term separately.

# Term 5: int cos(s)*r*f1^2*(f1')^2 ds and int sin(s)*r*f1^2*(f1')^2 ds
print("\n  --- Source term decomposition for f4 equation ---")

# Term 5: int_0^inf cos(s)*s*f1(s)^2*(f1'(s))^2 ds
# Note: the source is r * [-...], so the sign needs tracking.
# src4 = r * [-(term1 + term2 + term3 + term4 + term5)]
# where term5 = f1^2*(f1')^2 (positive sign in original, negative in source)

# So: contribution to B4 from term5 = -int cos(s)*s*f1^2*(f1')^2 ds

def term5(s):
    return s * f1(s)**2 * f1p(s)**2

T5c = mp.quadosc(lambda s: mp.cos(s)*term5(s), [0, mp.inf], omega=1)
T5s = mp.quadosc(lambda s: mp.sin(s)*term5(s), [0, mp.inf], omega=1)
print(f"  T5c = int cos*s*f1^2*(f1')^2 = {mp.nstr(T5c, 25)}")
print(f"  T5s = int sin*s*f1^2*(f1')^2 = {mp.nstr(T5s, 25)}")

# Term 3: -2*f1*f1'*f2'  → contribution: +2*int cos*s*f1*f1'*f2' ds (double negative)
# This is similar to K_c^(II) but with an extra f1 factor.
def term3_integrand_c(s):
    return mp.cos(s) * s * f1(s) * f1p(s) * f2p(s)
def term3_integrand_s(s):
    return mp.sin(s) * s * f1(s) * f1p(s) * f2p(s)

print("\n  Computing term3 integrals (f1*f1'*f2')...")
T3c = mp.quad(term3_integrand_c, [0, 50], maxdegree=8)
T3s = mp.quad(term3_integrand_s, [0, 50], maxdegree=8)
print(f"  T3c (R=50) = {mp.nstr(T3c, 20)}")
print(f"  T3s (R=50) = {mp.nstr(T3s, 20)}")

# Term 4: -f2*(f1')^2  → contribution: +int cos*s*f2*(f1')^2 ds
def f2_val(s):
    if s < mp.mpf('1e-10'): return mp.mpf(0)
    u2, _ = u2_and_u2p(s)
    return u2 / s

def term4_integrand_c(s):
    return mp.cos(s) * s * f2_val(s) * f1p(s)**2
def term4_integrand_s(s):
    return mp.sin(s) * s * f2_val(s) * f1p(s)**2

print("\n  Computing term4 integrals (f2*(f1')^2)...")
T4c = mp.quad(term4_integrand_c, [0, 50], maxdegree=8)
T4s = mp.quad(term4_integrand_s, [0, 50], maxdegree=8)
print(f"  T4c (R=50) = {mp.nstr(T4c, 20)}")
print(f"  T4s (R=50) = {mp.nstr(T4s, 20)}")

# Term 1: (f2')^2  → contribution: -int cos*s*(f2')^2 ds
def term1_integrand_c(s):
    fp2 = f2p(s)
    return mp.cos(s) * s * fp2**2
def term1_integrand_s(s):
    fp2 = f2p(s)
    return mp.sin(s) * s * fp2**2

print("\n  Computing term1 integrals ((f2')^2)...")
T1c = mp.quad(term1_integrand_c, [0, 50], maxdegree=8)
T1s = mp.quad(term1_integrand_s, [0, 50], maxdegree=8)
print(f"  T1c (R=50) = {mp.nstr(T1c, 20)}")
print(f"  T1s (R=50) = {mp.nstr(T1s, 20)}")

# Term 2: 2*f1'*f3' → This is the HARD one (requires f3)
# For now, let me compute the other terms and see what's left.

# B4 = -(T1c + T2c - 2*T3c - T4c + T5c) where T2c involves f3
# A4 = (T1s + T2s - 2*T3s - T4s + T5s)  (negated for cos amp)

# From terms 1,3,4,5 only (excluding term 2):
B4_partial = -(T1c - 2*T3c - T4c + T5c)
A4_partial = (T1s - 2*T3s - T4s + T5s)
print(f"\n  B4 (partial, excl f3 term, R=50) = {mp.nstr(B4_partial, 15)}")
print(f"  A4 (partial, excl f3 term, R=50) = {mp.nstr(A4_partial, 15)}")

# --- Summary ---
print(f"\n{'='*78}")
print(f"  SUMMARY")
print(f"{'='*78}")
print(f"  PROVEN: K_s^(I) = pi/12")
print(f"  K_s^(II) = {mp.nstr(KsII_full, 20)}")
print(f"  P_sin = pi/12 - 2*K_s^(II) = {mp.nstr(Psin, 20)}")
print(f"  A3 = -P_sin = {mp.nstr(A3_val, 20)}")
print(f"\n  For alpha_5 formula:")
print(f"    Known:  B2 = {mp.nstr(B2_val := mp.mpf(1)/2 - ln3/8, 20)}")
print(f"    Known:  A2 = {mp.nstr(A2_val := pi/8, 20)}")
print(f"    Known:  B3 = {mp.nstr(B3_val := mp.mpf('0.012615939290114711837850'), 20)}")
print(f"    NEW:    A3 = {mp.nstr(A3_val, 20)}")
print(f"    NEEDED: B4 (Q_cos) — requires f3 via double swap")
print(f"    NEEDED: A4 (-Q_sin) — requires f3 via double swap")
print(f"\n  Partial alpha_5 (without f3-dependent terms):")
# alpha_5 = q^2/2 + p*r - a*p*q + p^2*(a^2-b)/2 - p^4/8
# = A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
a = float(B2_val); b = float(B3_val); p = float(A2_val); q = float(A3_val)
a5_partial = q**2/2 - a*p*q + p**2*(a**2-b)/2 - p**4/8
print(f"    alpha_5 (excl A2*A4 term) = {a5_partial:.10f}")
print(f"    Missing: A2*A4 term contributes ~{p*float(A4_partial):.6f}")
print("="*78)
