#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c7_KcII_high_precision.py: High-precision numerical verification of
K_c^(II) = (ln 2 - 1)/12 - 9 pi^2/14080 via the swap decomposition.

Strategy: dps=40, quadosc with forced large period count.
Also try splitting integration range into [0,L] + [L,inf].
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
mp.mp.dps = 40

pi = mp.pi
ln2 = mp.log(2)
Ci = mp.ci
Si = mp.si

def f_fun(s):
    if s < mp.mpf('1e-20'): return mp.mpf(1) - s**2/6
    return mp.sin(s)/s

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*mp.cos(s) - mp.sin(s))/s**2

def A_kernel(t):
    return t * fp(t)**2

def Phi1(t):
    return -Ci(t)/2 + Ci(3*t)/2 - (mp.sin(t) + mp.sin(3*t))/(4*t)

def Phi2(t):
    return (Si(3*t) - Si(t))/2 - (mp.cos(t) - mp.cos(3*t))/(4*t)

def Phi3(t):
    return (3*mp.sin(t) - mp.sin(3*t))/(8*t) + 3*(Ci(3*t) - Ci(t))/8 + (mp.cos(3*t) - mp.cos(t))/(8*t**2)

def Phi4(t):
    return (5*mp.cos(t) - mp.cos(3*t))/(8*t) - pi/8 + 5*Si(t)/8 - 3*Si(3*t)/8 - (mp.sin(t) + mp.sin(3*t))/(8*t**2)

def A1_int(t):
    return mp.cos(t) * A_kernel(t) * Phi1(t)
def A2_int(t):
    return mp.sin(t) * A_kernel(t) * Phi2(t)
def A3_int(t):
    return mp.cos(t) * A_kernel(t) * Phi3(t)
def A4_int(t):
    return mp.sin(t) * A_kernel(t) * Phi4(t)

print("="*72)
print(f"  K_c^(II) HIGH PRECISION (dps={mp.mp.dps})")
print("="*72)

# Target values
KcI_th = (ln2 - 1)/6
KcII_target = (ln2 - 1)/12 - 9*pi**2/14080
Pcos_target = 9*pi**2/7040

print(f"\n  Predicted:")
print(f"    K_c^(I)     = (ln 2 - 1)/6 = {mp.nstr(KcI_th, 30)}")
print(f"    K_c^(II)    = (ln 2 - 1)/12 - 9 pi^2/14080")
print(f"                = {mp.nstr(KcII_target, 30)}")
print(f"    P_cos       = 9 pi^2/7040  = {mp.nstr(Pcos_target, 30)}")

print(f"\n  Computing A_1..A_4 via quadosc (dps={mp.mp.dps}, period=2*pi)...")
print(f"  This may take a minute...")

A1 = mp.quadosc(A1_int, [0, mp.inf], period=2*pi)
print(f"    A_1 = {mp.nstr(A1, 30)}")

A2 = mp.quadosc(A2_int, [0, mp.inf], period=2*pi)
print(f"    A_2 = {mp.nstr(A2, 30)}")

A3 = mp.quadosc(A3_int, [0, mp.inf], period=2*pi)
print(f"    A_3 = {mp.nstr(A3, 30)}")

A4 = mp.quadosc(A4_int, [0, mp.inf], period=2*pi)
print(f"    A_4 = {mp.nstr(A4, 30)}")

KcII_num = -A1 - A2 + A3 - A4
diff_KcII = KcII_num - KcII_target

print(f"\n  RESULT:")
print(f"    K_c^(II) num    = {mp.nstr(KcII_num, 30)}")
print(f"    K_c^(II) target = {mp.nstr(KcII_target, 30)}")
print(f"    diff            = {mp.nstr(diff_KcII, 5)}")
print(f"    agreement: {abs(float(diff_KcII)):.2e}")

Pcos_num = KcI_th - 2*KcII_num
diff_Pcos = Pcos_num - Pcos_target
print(f"\n    P_cos num    = {mp.nstr(Pcos_num, 30)}")
print(f"    P_cos target = {mp.nstr(Pcos_target, 30)}")
print(f"    diff         = {mp.nstr(diff_Pcos, 5)}")
print(f"    agreement: {abs(float(diff_Pcos)):.2e}")

print(f"\n{'='*72}")
