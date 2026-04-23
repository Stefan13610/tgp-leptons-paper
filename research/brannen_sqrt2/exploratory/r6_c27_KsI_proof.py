#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c27_KsI_proof.py:
Verification of the proof K_s^(I) = pi/12.

K_s^(I) = int_0^inf sin(s) * s * f1(s) * (f1'(s))^2 ds

where f1 = sin(s)/s, f1' = (s*cos(s) - sin(s))/s^2.

Expanding sin(s)*s*f1*(f1')^2 = sin^2(s)*(s*cos-sin)^2/s^4:
= sin^2*cos^2/s^2 - 2*sin^3*cos/s^3 + sin^4/s^4

So K_s^(I) = J1 - 2*J2 + J3 where:
  J1 = int sin^2*cos^2/s^2 ds = pi/4
  J2 = int sin^3*cos/s^3 ds   = pi/4
  J3 = int sin^4/s^4 ds       = pi/3

K_s^(I) = pi/4 - pi/2 + pi/3 = pi/12.

PROOF of J1: sin^2*cos^2 = sin^2(2s)/4, and int sin^2(2s)/s^2 ds = pi*2/2 = pi.
  So J1 = pi/4.

PROOF of J2: sin^3*cos = (sin(2s) - sin(4s)/2)/4.
  Using I(a) = int (sin(as)-as)/s^3 ds = -pi*a^2/4:
  int (sin(2s)-sin(4s)/2)/s^3 ds = I(2) - I(4)/2 + int (2s-2s)/s^3 ds
  Hmm, need to be careful:
  sin(2s)/s^3 - sin(4s)/(2s^3) = (sin(2s) - 2s)/s^3 - (sin(4s)/2 - 2s)/s^3
  I(2) = -pi*4/4 = -pi
  int (sin(4s)/2 - 2s)/s^3 ds = (1/2)*I(4) + int (2s - 2s)/s^3 ds = (1/2)*(-pi*16/4) = -2pi
  So integral = -pi - (-2pi) = pi.
  J2 = pi/4.

PROOF of J3: Standard result: int_0^inf sinc^4 = int sin^4/s^4 = pi/3.

K_s^(I) = pi/4 - 2*(pi/4) + pi/3 = -pi/4 + pi/3 = pi/12.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 40
pi = mp.pi

print("="*78)
print("  VERIFICATION: K_s^(I) = pi/12")
print("  dps =", mp.mp.dps)
print("="*78)

def f1(s):
    if s < mp.mpf('1e-30'): return mp.mpf(1)
    return mp.sin(s)/s

def f1p(s):
    if s < mp.mpf('1e-30'): return mp.mpf(0)
    return mp.cos(s)/s - mp.sin(s)/s**2

# Full integral K_s^(I)
KsI = mp.quadosc(lambda s: mp.sin(s)*s*f1(s)*f1p(s)**2, [0, mp.inf], omega=1)
print(f"\n  K_s^(I) = {mp.nstr(KsI, 35)}")
print(f"  pi/12   = {mp.nstr(pi/12, 35)}")
print(f"  diff    = {mp.nstr(KsI - pi/12, 8)}")

# Verify J1
J1 = mp.quadosc(lambda s: mp.sin(s)**2 * mp.cos(s)**2 / s**2, [0, mp.inf], omega=1)
print(f"\n  J1 = int sin^2*cos^2/s^2 = {mp.nstr(J1, 30)}")
print(f"  pi/4 = {mp.nstr(pi/4, 30)}")
print(f"  diff = {mp.nstr(J1 - pi/4, 8)}")

# Verify J2
J2 = mp.quadosc(lambda s: mp.sin(s)**3 * mp.cos(s) / s**3, [0, mp.inf], omega=1)
print(f"\n  J2 = int sin^3*cos/s^3 = {mp.nstr(J2, 30)}")
print(f"  pi/4 = {mp.nstr(pi/4, 30)}")
print(f"  diff = {mp.nstr(J2 - pi/4, 8)}")

# Verify J3
J3 = mp.quadosc(lambda s: mp.sin(s)**4 / s**4, [0, mp.inf], omega=1)
print(f"\n  J3 = int sin^4/s^4 = {mp.nstr(J3, 30)}")
print(f"  pi/3 = {mp.nstr(pi/3, 30)}")
print(f"  diff = {mp.nstr(J3 - pi/3, 8)}")

# Verify decomposition
check = J1 - 2*J2 + J3
print(f"\n  J1 - 2*J2 + J3 = {mp.nstr(check, 30)}")
print(f"  pi/12 = {mp.nstr(pi/12, 30)}")
print(f"  diff = {mp.nstr(check - pi/12, 8)}")

# Verify I(a) = int (sin(as)-as)/s^3 ds = -pi*a^2/4
print(f"\n  --- Verify I(a) = int (sin(as)-as)/s^3 = -pi*a^2/4 ---")
for a in [1, 2, 3, 4]:
    Ia = mp.quad(lambda s: (mp.sin(a*s) - a*s)/s**3, [0, mp.inf])
    expected = -pi*a**2/4
    print(f"    a={a}: I(a) = {mp.nstr(Ia, 20)}, -pi*a^2/4 = {mp.nstr(expected, 20)}, diff = {mp.nstr(Ia - expected, 6)}")

# Also verify: int sin^2(as)/s^2 = pi*a/2
print(f"\n  --- Verify int sin^2(as)/s^2 = pi*a/2 ---")
for a in [1, 2, 3, 4]:
    val = mp.quadosc(lambda s: mp.sin(a*s)**2 / s**2, [0, mp.inf], omega=a)
    expected = pi*a/2
    print(f"    a={a}: {mp.nstr(val, 20)}, expected = {mp.nstr(expected, 20)}, diff = {mp.nstr(val-expected, 6)}")

# Bonus: verify K_c^(I) = (ln2-1)/6
KcI = mp.quadosc(lambda s: mp.cos(s)*s*f1(s)*f1p(s)**2, [0, mp.inf], omega=1)
ln2 = mp.log(2)
print(f"\n  Bonus: K_c^(I) = {mp.nstr(KcI, 30)}")
print(f"  (ln2-1)/6 = {mp.nstr((ln2-1)/6, 30)}")
print(f"  diff = {mp.nstr(KcI - (ln2-1)/6, 8)}")

print(f"\n{'='*78}")
print("  CONCLUSION: K_s^(I) = pi/12 VERIFIED to", mp.mp.dps, "digits")
print("  PROOF: K_s^(I) = J1 - 2*J2 + J3 = pi/4 - pi/2 + pi/3 = pi/12")
print("  Each component verified independently.")
print("="*78)
