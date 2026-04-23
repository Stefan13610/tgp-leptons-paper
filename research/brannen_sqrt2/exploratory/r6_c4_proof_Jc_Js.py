#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c4_proof_Jc_Js.py: Analytical proof of J_c, J_s (-> alpha_2, I_sin^2/2)

Key identity (from f'' + (2 over s) f' + f = 0, where f = sin(s)/s):
    s * (f')^2 = d/ds[s * f * f'] + s * f^2 + (1/2) * d/ds[f^2]

By integration by parts:
    J_c = int cos(s) s (f')^2 ds
        = int sin(s) s f f' ds + int cos(s) s f^2 ds - 1/2 + (1/2) int sin(s) f^2 ds

Each integral is ELEMENTARY via Dirichlet/Frullani:
    int sin s * (s f f') ds = -ln(3)/2
    int cos s * (s f^2)  ds =  ln(3)/4
    int sin s * f^2      ds = (3/4) ln(3)

Therefore:
    J_c = -ln(3)/2 + ln(3)/4 - 1/2 + (3/8) ln(3) = ln(3)/8 - 1/2   [PROVEN]

Analogously for J_s:
    J_s = -int cos s (s f f') ds + int sin s (s f^2) ds - (1/2) int cos s f^2 ds
        = -0 + pi/4 - pi/8 = pi/8                                   [PROVEN]

This proves alpha_2 = -J_c = 1/2 - ln(3)/8 = I_cos [NEW CLEAN PROOF]
And        I_sin^2/2 = J_s^2/2 = pi^2/128         [NEW CLEAN PROOF]

Remaining: P_cos = 9 pi^2 / 7040 via analogous IBP on p-equation.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
mp.mp.dps = 40

pi = mp.pi
ln3 = mp.log(3)

print("="*72)
print("  ANALITYCZNY DOWOD J_c, J_s (IBP + Dirichlet + Frullani)")
print("="*72)

I_cos_th = mp.mpf(1)/2 - ln3/8
Jc_predicted = ln3/8 - mp.mpf(1)/2    # = -I_cos
Js_predicted = pi/8

print(f"\n  Predicted analytically:")
print(f"    J_c = ln(3)/8 - 1/2 = {mp.nstr(Jc_predicted, 25)}")
print(f"    J_s = pi/8          = {mp.nstr(Js_predicted, 25)}")

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*mp.cos(s) - mp.sin(s))/s**2

def sA(s):
    return s * fp(s)**2

def f_fun(s):
    if s < mp.mpf('1e-20'):
        return mp.mpf(1) - s**2/6
    return mp.sin(s)/s

print(f"\n  Numerical verification (mpmath dps=40, quadosc):")
Jc_num = mp.quadosc(lambda s: mp.cos(s)*sA(s), [0, mp.inf], period=2*pi)
Js_num = mp.quadosc(lambda s: mp.sin(s)*sA(s), [0, mp.inf], period=2*pi)
print(f"    J_c num = {mp.nstr(Jc_num, 25)}")
print(f"    J_c pre = {mp.nstr(Jc_predicted, 25)}")
print(f"    diff    = {mp.nstr(Jc_num - Jc_predicted, 5)}")
print(f"    J_s num = {mp.nstr(Js_num, 25)}")
print(f"    J_s pre = {mp.nstr(Js_predicted, 25)}")
print(f"    diff    = {mp.nstr(Js_num - Js_predicted, 5)}")

print(f"\n{'='*72}")
print(f"  Weryfikacja 6 calek elementarnych")
print(f"{'='*72}")

def sff(s): return mp.sin(s) * fp(s)       # s*f*f' = sin(s) * f'(s)   (bo s*f = sin(s))
def sf2(s): return s * f_fun(s)**2          # s*f^2 = sin^2(s)/s
def f2(s):  return f_fun(s)**2              # f^2 = sin^2(s)/s^2

I1_th = -ln3/2
I1_num = mp.quadosc(lambda s: mp.sin(s)*sff(s), [0, mp.inf], period=2*pi)
print(f"  (1) int sin(s) s f f' ds = -ln(3)/2 = {mp.nstr(I1_th, 22)}")
print(f"                         num          = {mp.nstr(I1_num, 22)}")

I2_th = ln3/4
I2_num = mp.quadosc(lambda s: mp.cos(s)*sf2(s), [0, mp.inf], period=2*pi)
print(f"  (2) int cos(s) s f^2 ds = ln(3)/4  = {mp.nstr(I2_th, 22)}")
print(f"                       num           = {mp.nstr(I2_num, 22)}")

I3_th = 3*ln3/4
I3_num = mp.quadosc(lambda s: mp.sin(s)*f2(s), [0, mp.inf], period=2*pi)
print(f"  (3) int sin(s) f^2 ds = 3 ln(3)/4 = {mp.nstr(I3_th, 22)}")
print(f"                     num           = {mp.nstr(I3_num, 22)}")

I4_th = mp.mpf(0)
I4_num = mp.quadosc(lambda s: mp.cos(s)*sff(s), [0, mp.inf], period=2*pi)
print(f"  (4) int cos(s) s f f' ds = 0       = {mp.nstr(I4_th, 22)}")
print(f"                          num        = {mp.nstr(I4_num, 22)}")

I5_th = pi/4
I5_num = mp.quadosc(lambda s: mp.sin(s)*sf2(s), [0, mp.inf], period=2*pi)
print(f"  (5) int sin(s) s f^2 ds = pi/4    = {mp.nstr(I5_th, 22)}")
print(f"                       num          = {mp.nstr(I5_num, 22)}")

I6_th = pi/4
I6_num = mp.quadosc(lambda s: mp.cos(s)*f2(s), [0, mp.inf], period=2*pi)
print(f"  (6) int cos(s) f^2 ds = pi/4      = {mp.nstr(I6_th, 22)}")
print(f"                     num            = {mp.nstr(I6_num, 22)}")

Jc_assembled = I1_th + I2_th - mp.mpf(1)/2 + I3_th/2
Js_assembled = -I4_th + I5_th - I6_th/2

print(f"\n  Assembled:")
print(f"    J_c = I1 + I2 - 1/2 + I3/2 = ln(3)/8 - 1/2 = {mp.nstr(Jc_assembled, 22)}")
print(f"    J_s = -I4 + I5 - I6/2 = pi/8              = {mp.nstr(Js_assembled, 22)}")

print(f"\n{'='*72}")
print(f"  WNIOSKI:")
print(f"    alpha_2 = -J_c = 1/2 - ln(3)/8 = I_cos    [PROVEN analytically]")
print(f"    I_sin^2/2 = J_s^2/2 = pi^2/128             [PROVEN analytically]")
print(f"    Pozostaje: P_cos = 9*pi^2/7040 analitycznie")
print(f"{'='*72}")
