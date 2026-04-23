#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c11_alpha3_ultra_precision.py: Ultra-high-precision numerical check of
K_c^(II) and alpha_3 by computing P_cos via the ORIGINAL formula (not through h).

P_cos = int_0^inf cos(r) r [f(r)(f'(r))^2 - 2 f'(r) h'(r)] dr

Strategy:
1. Solve h-equation: h'' + (2/r)h' + h = -(f')^2, h(0) = 0, h bounded near origin.
   Use Green's function directly: h(r) = -u(r)/r where u(r) = int_0^r sin(r-s) s (f')^2 ds.
   Write u(r) = sin(r) J_c(r) - cos(r) J_s(r).

2. Evaluate P_cos numerically to 12+ digits. Compare to 9 pi^2/7040.

Uses mpmath with dps=50 for inner evaluations. This is slow but converges.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
mp.mp.dps = 50

pi = mp.pi
ln2 = mp.log(2)

print("="*72)
print(f"  ULTRA-PRECISION CHECK: alpha_3 and Kc^(II)  (dps={mp.mp.dps})")
print("="*72)

def f_fun(s):
    if s < mp.mpf('1e-20'): return mp.mpf(1) - s**2/6 + s**4/120
    return mp.sin(s)/s

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360 - s**9/3991680
    return (s*mp.cos(s) - mp.sin(s))/s**2

Ci = mp.ci
Si = mp.si

# Phi_i closed forms (already verified to 30+ digits)
def Phi1(t):
    return -Ci(t)/2 + Ci(3*t)/2 - (mp.sin(t) + mp.sin(3*t))/(4*t)

def Phi2(t):
    return (Si(3*t) - Si(t))/2 - (mp.cos(t) - mp.cos(3*t))/(4*t)

def Phi3(t):
    return (3*mp.sin(t) - mp.sin(3*t))/(8*t) + 3*(Ci(3*t) - Ci(t))/8 + (mp.cos(3*t) - mp.cos(t))/(8*t**2)

def Phi4(t):
    return (5*mp.cos(t) - mp.cos(3*t))/(8*t) - pi/8 + 5*Si(t)/8 - 3*Si(3*t)/8 - (mp.sin(t) + mp.sin(3*t))/(8*t**2)

def A_ker(t): return t * fp(t)**2

def A1_int(t): return mp.cos(t) * A_ker(t) * Phi1(t)
def A2_int(t): return mp.sin(t) * A_ker(t) * Phi2(t)
def A3_int(t): return mp.cos(t) * A_ker(t) * Phi3(t)
def A4_int(t): return mp.sin(t) * A_ker(t) * Phi4(t)

# Compute with quadosc (may take a while at dps=50)
print(f"\n  Computing A_1, A_2, A_3, A_4 at dps={mp.mp.dps}...")

A1 = mp.quadosc(A1_int, [0, mp.inf], period=2*pi)
print(f"    A_1 = {mp.nstr(A1, 40)}")

A2 = mp.quadosc(A2_int, [0, mp.inf], period=2*pi)
print(f"    A_2 = {mp.nstr(A2, 40)}")

A3 = mp.quadosc(A3_int, [0, mp.inf], period=2*pi)
print(f"    A_3 = {mp.nstr(A3, 40)}")

A4 = mp.quadosc(A4_int, [0, mp.inf], period=2*pi)
print(f"    A_4 = {mp.nstr(A4, 40)}")

KcII = -A1 - A2 + A3 - A4
KcI = (ln2 - 1)/6
Pcos = KcI - 2*KcII
alpha3 = pi**2/128 + Pcos

KcII_target = (ln2 - 1)/12 - 9*pi**2/14080
Pcos_target = 9*pi**2/7040
alpha3_target = pi**2/110

print(f"\n{'='*72}")
print(f"  RESULTS (ultra-precision):")
print(f"{'='*72}")
print(f"    K_c^(II)   = {mp.nstr(KcII, 35)}")
print(f"    target     = {mp.nstr(KcII_target, 35)}")
print(f"    diff       = {mp.nstr(KcII - KcII_target, 8)}")
print(f"    P_cos      = {mp.nstr(Pcos, 35)}")
print(f"    target     = {mp.nstr(Pcos_target, 35)}")
print(f"    diff       = {mp.nstr(Pcos - Pcos_target, 8)}")
print(f"    alpha_3    = {mp.nstr(alpha3, 35)}")
print(f"    target     = {mp.nstr(alpha3_target, 35)}")
print(f"    diff       = {mp.nstr(alpha3 - alpha3_target, 8)}")

# Look for alternative closed forms of alpha_3:
# Candidates: pi^2/110, (some simple rational)*pi^2, possibly with ln terms.
# Check pi^2/(integer) for integer ~ 110
print(f"\n{'='*72}")
print(f"  CANDIDATE MATCH CHECK for alpha_3:")
print(f"{'='*72}")
for n in range(108, 115):
    val = pi**2/n
    print(f"    pi^2/{n} = {mp.nstr(val, 35)}   diff = {mp.nstr(alpha3 - val, 8)}")

print("="*72)
