#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c8_moment_integrals.py: Table of moment integrals needed for A_1, A_2, A_3, A_4.

Each A_i expands (via frequency decomposition of trig * (f')^2 and Phi_i) into a
sum of elementary integrals like:
    M1(a,b) = int_0^inf sin(a t) Si(b t)/t dt
    M2(a,b) = int_0^inf sin(a t) Ci(b t)/t dt
    M3(a,b) = int_0^inf cos(a t) Si(b t)/t dt
    M4(a,b) = int_0^inf cos(a t) Ci(b t)/t dt
    plus powers 1/t^2, 1/t^3, 1/t^4 (reducible by IBP).

Known identities:
    M1(a,b) + M1(b,a) = pi^2/4        (by IBP sym)
    M1(a,a) = pi^2/8
    M1(a,b) for a > b: = chi_2(b/a) = (1/2)[Li_2(b/a) - Li_2(-b/a)]
                      where chi_2(x) = x + x^3/9 + x^5/25 + ...

This file:
 - Verifies numerically a set of reference values
 - Provides the chi_2(1/3) and similar constants
 - Sketches the derivation needed for K_c^(II).
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
mp.mp.dps = 30

pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)
Ci = mp.ci
Si = mp.si

print("="*72)
print("  MOMENT INTEGRAL TABLE")
print("="*72)

# M1(a, b) = int_0^inf sin(a t) Si(b t)/t dt
# For a > b > 0: M1(a,b) = chi_2(b/a) = sum (b/a)^(2k+1)/(2k+1)^2
# For a = b:     M1(a,a) = pi^2/8
# For a < b:     M1(a,b) = pi^2/4 - M1(b,a)  [by IBP symmetry]

def chi2(x):
    """Legendre chi function chi_2(x) = sum x^(2k+1)/(2k+1)^2"""
    return (mp.polylog(2, x) - mp.polylog(2, -x))/2

def M1(a, b):
    """M1(a,b) = int_0^inf sin(a t) Si(b t)/t dt"""
    if a == b:
        return pi**2/8
    if a > b:
        return chi2(mp.mpf(b)/a)
    else:  # a < b
        return pi**2/4 - chi2(mp.mpf(a)/b)

# ==== Verify M1(a,b) numerically ====
print("\n  M1(a,b) = int sin(a*t) Si(b*t)/t dt:")
for a, b in [(1, 1), (1, 3), (3, 1), (3, 3), (1, 2), (2, 1)]:
    th = M1(a, b)
    num = mp.quadosc(lambda t: mp.sin(a*t) * Si(b*t) / t, [mp.mpf('1e-30'), mp.inf], period=2*pi/max(a,b))
    print(f"    M1({a},{b}) = {mp.nstr(th, 18):>22}  num={mp.nstr(num, 18):>22}  diff={mp.nstr(th-num, 4)}")

# ==== M2(a, b) = int sin(a t) Ci(b t)/t dt ====
# IBP: u = Ci(bt), du = cos(bt)/t dt; dv = sin(at)/t dt, v = Si(at)
# M2(a,b) = [Si(at) Ci(bt)]_0^inf - int Si(at) cos(bt)/t dt
# At infty: Si -> pi/2, Ci -> 0. Product: 0.
# At 0: Si(at)=at + O, Ci(bt) = gamma + ln(bt) + O. Si*Ci = at(gamma + ln bt) -> 0.
# So M2(a,b) = -int Si(at) cos(bt)/t dt.

# Also: M2(a,b) = -int Si(at) cos(bt)/t dt. And by parts another way:
# Use: d/dt[Si(at) Ci(bt)] = a cos(at) Ci(bt)/... hmm wait d Si/dt = sin(at)/t... no:
# Si(x) = int_0^x sin(u)/u du, so d/dt Si(at) = sin(at) · (a/at) = sin(at)/t. Good.
# d/dt Ci(bt) = cos(bt)/t.

# Known: int_0^inf Si(at) cos(bt)/t dt: for a,b > 0:
# = pi/4 * ln|(a+b)/(a-b)| + ???
# Standard result (see Gradshteyn 6.233 or similar): need to look up.

# Try numerical M2(a,b):
print("\n  M2(a,b) = int sin(a*t) Ci(b*t)/t dt:")
for a, b in [(1, 1), (1, 3), (3, 1)]:
    num = mp.quadosc(lambda t: mp.sin(a*t) * Ci(b*t) / t, [mp.mpf('1e-30'), mp.inf], period=2*pi/max(a,b))
    print(f"    M2({a},{b}) num = {mp.nstr(num, 20)}")

# ==== Verify key constant chi_2(1/3) ====
c13 = chi2(mp.mpf(1)/3)
print(f"\n  chi_2(1/3) = {mp.nstr(c13, 25)}")
print(f"    check: sum = 1/3 + 1/(9*27) + 1/(25*243) + ...")
s = mp.mpf(0)
x = mp.mpf(1)/3
for k in range(0, 50):
    s += x**(2*k+1)/(2*k+1)**2
print(f"    partial sum (50 terms) = {mp.nstr(s, 25)}")

# ==== Catalan-type constants ====
print(f"\n  chi_2(1) = pi^2/8 = {mp.nstr(pi**2/8, 25)}")
print(f"  chi_2(1/2) = {mp.nstr(chi2(mp.mpf(1)/2), 25)}")
print(f"  chi_2(1/3) = {mp.nstr(chi2(mp.mpf(1)/3), 25)}")

# ==== Useful: pi^2/110 vs pi^2/128 and the gap ====
Pcos_target = 9*pi**2/7040
print(f"\n  P_cos target = 9 pi^2/7040 = {mp.nstr(Pcos_target, 25)}")
print(f"  alpha_3 target = pi^2/110  = {mp.nstr(pi**2/110, 25)}")

# KcII predicted
KcII_target = (ln2 - 1)/12 - 9*pi**2/14080
print(f"  K_c^(II) predicted = {mp.nstr(KcII_target, 25)}")

# Does K_c^(II) decompose nicely? A possible clean form:
# K_c^(II) should equal some combination of chi_2(1/3), pi^2/k, ln3, ln2.
#
# Check: (ln 2 - 1)/12 - 9 pi^2/14080
# 14080 = 128 * 110 = 2^7 * 5 * 11 * ... = 2^7 * 110
# Actually 14080 = 128 * 110 = 14080. And 110 = 2*5*11. So 14080 = 2^8 * 5 * 11. Yes.

print(f"\n  Factor check: 14080 = {2**7 * 110} = 2^7 * 110 = 2^8 * 5 * 11")
print(f"               7040  = {2**6 * 110} = 2^6 * 110 = 2^7 * 5 * 11")
print(f"               110   = 2*5*11")

print("="*72)
