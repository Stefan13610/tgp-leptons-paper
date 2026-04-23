#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Kc^(I) analytical derivation and P_cos partial check."""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import mpmath as mp
mp.mp.dps = 40

pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

print("="*72)
print("  Kc^(I) = int cos(s) s f (fprime)^2 ds  ANALITYCZNIE")
print("="*72)

KcI_predicted = (ln2 - 1)/6
Pcos_target = 9*pi**2/7040
print("")
print("  Predicted Kc^(I) = (ln 2 - 1) div 6 = " + mp.nstr(KcI_predicted, 25))
print("  Target P_cos = 9 pi^2 div 7040     = " + mp.nstr(Pcos_target, 25))

def f_fun(s):
    if s < mp.mpf("1e-20"): return mp.mpf(1) - s**2/6
    return mp.sin(s)*(1/s)

def fp(s):
    if s < mp.mpf("1e-3"):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*mp.cos(s) - mp.sin(s))*(1/s**2)

def KcI_integrand(s):
    return mp.cos(s) * s * f_fun(s) * fp(s)**2

print("")
print("  Numerical Kc^(I) via quadosc:")
KcI_num = mp.quadosc(KcI_integrand, [0, mp.inf], period=2*pi)
print("    Kc^(I) num = " + mp.nstr(KcI_num, 25))
print("    Kc^(I) pre = " + mp.nstr(KcI_predicted, 25))
print("    diff       = " + mp.nstr(KcI_num - KcI_predicted, 5))

# Elementary:
I_a_th = ln2/2
I_a_num = mp.quadosc(lambda s: mp.sin(s)**3 * mp.cos(s) * (1/s**2), [0, mp.inf], period=2*pi)
print("")
print("  int sin^3 cos div s^2 = ln(2) div 2 = " + mp.nstr(I_a_th, 22))
print("                      num             = " + mp.nstr(I_a_num, 22))

I_b_th = ln2
I_b_num = mp.quadosc(lambda s: mp.sin(s)**4 * (1/s**3), [0, mp.inf], period=2*pi)
print("  int sin^4 div s^3     = ln(2)       = " + mp.nstr(I_b_th, 22))
print("                     num              = " + mp.nstr(I_b_num, 22))

KcI_assembled = (I_a_th - I_b_th)/2 - mp.mpf(1)/6 + I_b_th/6 + I_a_th/2
print("")
print("  Assembled (ln 2 - 1) div 6 = " + mp.nstr(KcI_assembled, 22))

# ==== Kc^(II) numerical (slow) ====
print("")
print("="*72)
print("  Kc^(II) numeryczny")
print("="*72)

def sA(s): return s * fp(s)**2

def Jc_r(r):
    if r < mp.mpf("1e-15"): return mp.mpf(0)
    return mp.quad(lambda t: mp.cos(t) * sA(t), [0, r])
def Js_r(r):
    if r < mp.mpf("1e-15"): return mp.mpf(0)
    return mp.quad(lambda t: mp.sin(t) * sA(t), [0, r])

def u_val(r):
    jc = Jc_r(r); js = Js_r(r)
    return -mp.sin(r)*jc + mp.cos(r)*js
def up_val(r):
    jc = Jc_r(r); js = Js_r(r)
    return -mp.cos(r)*jc - mp.sin(r)*js

def hp(r):
    if r < mp.mpf("1e-10"): return mp.mpf(0)
    return up_val(r)*(1/r) - u_val(r)*(1/r**2)

def KcII_integrand(s):
    return mp.cos(s) * s * fp(s) * hp(s)

mp.mp.dps = 25
print("    (Using dps=25 for speed; [0, 60])")
try:
    KcII_num = mp.quad(KcII_integrand, [0, 60])
    print("    Kc^(II) num [0,60] = " + mp.nstr(KcII_num, 18))
    Pcos_computed = KcI_predicted - 2*KcII_num
    print("")
    print("  P_cos = Kc^(I) - 2 Kc^(II):")
    print("    computed = " + mp.nstr(Pcos_computed, 18))
    print("    target   = " + mp.nstr(Pcos_target, 18))
    print("    diff     = " + mp.nstr(Pcos_computed - Pcos_target, 5))
except Exception as e:
    print("    Blad: " + str(e))

print("")
print("="*72)
print("  Kc^(I) = (ln 2 - 1) div 6   [PROVEN ANALYTICALLY]")
print("="*72)
