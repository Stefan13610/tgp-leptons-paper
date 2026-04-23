#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c14_alpha3_precision_convergence.py:

Decisive precision-convergence test: run the swap-based alpha_3 at multiple dps
(25, 40, 60) and check whether the residual (alpha_3 - pi^2/110) stabilizes
or shrinks with higher precision.

If stable => alpha_3 is NOT pi^2/110 exactly; the conjecture is wrong.
If shrinking => quadosc convergence issue, real value is pi^2/110.

Also compute symbolic decompositions for A_1 as A_1^pure (no Ci) + A_1^Ci
separately, so we can tell which piece is the error source.
"""
import sys, io, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360 - s**9/3991680
    return (s*mp.cos(s) - mp.sin(s))/s**2

def run_at_dps(dps):
    mp.mp.dps = dps
    pi = mp.pi
    ln2 = mp.log(2)
    ln3 = mp.log(3)
    Ci = mp.ci
    Si = mp.si

    def Aker(t): return t * fp(t)**2
    def Phi1(t):
        return -Ci(t)/2 + Ci(3*t)/2 - (mp.sin(t) + mp.sin(3*t))/(4*t)
    def Phi2(t):
        return (Si(3*t) - Si(t))/2 - (mp.cos(t) - mp.cos(3*t))/(4*t)
    def Phi3(t):
        return (3*mp.sin(t) - mp.sin(3*t))/(8*t) + 3*(Ci(3*t) - Ci(t))/8 + (mp.cos(3*t) - mp.cos(t))/(8*t**2)
    def Phi4(t):
        return (5*mp.cos(t) - mp.cos(3*t))/(8*t) - pi/8 + 5*Si(t)/8 - 3*Si(3*t)/8 - (mp.sin(t) + mp.sin(3*t))/(8*t**2)

    print(f"\n  --- dps = {dps} ---")
    t0 = time.time()
    A1 = mp.quadosc(lambda t: mp.cos(t) * Aker(t) * Phi1(t), [0, mp.inf], period=2*pi)
    A2 = mp.quadosc(lambda t: mp.sin(t) * Aker(t) * Phi2(t), [0, mp.inf], period=2*pi)
    A3 = mp.quadosc(lambda t: mp.cos(t) * Aker(t) * Phi3(t), [0, mp.inf], period=2*pi)
    A4 = mp.quadosc(lambda t: mp.sin(t) * Aker(t) * Phi4(t), [0, mp.inf], period=2*pi)
    elapsed = time.time() - t0

    KcI = (ln2 - 1)/6
    KcII = -A1 - A2 + A3 - A4
    Pcos = KcI - 2*KcII
    alpha3 = pi**2/128 + Pcos
    target_alpha3 = pi**2/110
    diff = alpha3 - target_alpha3

    print(f"    A_1 = {mp.nstr(A1, min(30, dps-3))}")
    print(f"    A_2 = {mp.nstr(A2, min(30, dps-3))}")
    print(f"    A_3 = {mp.nstr(A3, min(30, dps-3))}")
    print(f"    A_4 = {mp.nstr(A4, min(30, dps-3))}")
    print(f"    Kc^(II) = {mp.nstr(KcII, min(30, dps-3))}")
    print(f"    alpha_3 = {mp.nstr(alpha3, min(30, dps-3))}")
    print(f"    diff (alpha_3 - pi^2/110) = {mp.nstr(diff, 12)}")
    print(f"    (took {elapsed:.1f} s)")
    return alpha3, diff

print("="*72)
print("  PRECISION-CONVERGENCE TEST for alpha_3 via swap decomposition")
print("="*72)

results = []
for dps in [25, 40, 60]:
    a, d = run_at_dps(dps)
    results.append((dps, a, d))

print(f"\n{'='*72}")
print("  SUMMARY:")
print(f"{'='*72}")
print(f"  {'dps':>6}  {'alpha_3':>36}  {'diff from pi^2/110':>22}")
for dps, a, d in results:
    print(f"  {dps:>6}  {mp.nstr(a, 32):>36}  {mp.nstr(d, 12):>22}")

# Verdict
diffs = [d for _, _, d in results]
if all(abs(float(d) + 1.45e-6) < 1e-7 for d in diffs):
    print("\n  VERDICT: diff stable at ~-1.45e-6 across all precisions.")
    print("          alpha_3 is NOT pi^2/110 exactly — conjecture revision needed.")
elif abs(float(diffs[-1])) < abs(float(diffs[0]))/10:
    print("\n  VERDICT: diff shrinking with precision — quadosc issue.")
    print("          alpha_3 = pi^2/110 likely correct; need higher-precision method.")
else:
    print("\n  VERDICT: diffs mixed — need deeper investigation.")

print("="*72)
