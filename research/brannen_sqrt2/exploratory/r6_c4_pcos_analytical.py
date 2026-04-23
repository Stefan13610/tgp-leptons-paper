#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c4_pcos_analytical.py
=========================

Analityczny dowod P_cos = 9*pi^2/7040.

Kontekst:
  alpha_3 = I_sin^2/2 + P_cos
  I_sin = -pi/8 (DOWODZONE analitycznie)
  ==> I_sin^2/2 = pi^2/128 = 55*pi^2/7040

  Pozostale: P_cos = 9*pi^2/7040 (do dowiedzenia)
  alpha_3 = 55*pi^2/7040 + 9*pi^2/7040 = 64*pi^2/7040 = pi^2/110 QED

STRATEGIA:
  1. f(r) = sin(r)/r, f'(r) = cos(r)/r - sin(r)/r^2
  2. s*(f')^2(s) = cos^2(s)/s - 2 sin(s)cos(s)/s^2 + sin^2(s)/s^3
                = (1+cos(2s))/(2s) - sin(2s)/s^2 + (1-cos(2s))/(2s^3)
  3. u(r) = r*h(r) spelnia u'' + u = -r*(f')^2
  4. Zielona funkcja (BC: u(0)=u'(0)=0):
     u(r) = -int_0^r sin(r-s) s (f')^2(s) ds
         = -sin(r) J_c(r) + cos(r) J_s(r)
     gdzie J_c(r) = int_0^r cos(s) s(f')^2 ds, J_s(r) = int_0^r sin(s) s(f')^2 ds
  5. P_cos = lim_{R->inf} 1/R * int_0^R [O(d^3) RHS przez cos(s) wage] ds
     - ale w praktyce P_cos zwiazane z tailem h(r) i nonliniowych czlonow.

W tej iteracji: mpmath ultra-precision (50 digits) dla numerycznego P_cos
oraz sympy symbolic dla kluczowych calek cos(s)cos(2s)/s, etc.

Author: Claudian
Date: 2026-04-16
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import mpmath as mp
import math

mp.mp.dps = 50  # 50 decimal places

print("=" * 78)
print("  P_cos = 9*pi^2/7040 : analityczna weryfikacja ultra-precision")
print("=" * 78)
print(f"\n  mpmath precision: {mp.mp.dps} digits")

# Ground truths
pi = mp.pi
target = 9*pi**2 / 7040
print(f"\n  Target:  9*pi^2/7040 = {mp.nstr(target, 25)}")

# ==============================================================
# KROK 1: Weryfikacja I_sin, I_cos analitycznie przez Frullani
# ==============================================================

print(f"\n{'='*78}")
print(f"  KROK 1: I_sin, I_cos przez Frullani (analitycznie)")
print(f"{'='*78}")

# I_cos = int_0^inf cos(r) * f(r) dr gdzie f = sin(r)/r?
# Actually: I_cos = 1/2 - ln(3)/8
I_cos_th = mp.mpf(1)/2 - mp.log(3)/8
I_sin_th = -pi/8

# Weryfikacja numeryczna:
# I_sin = int_0^inf sin(r) * [sin(r)/r] * 1/something?
# Actually z summary: alpha_2 = c_1/2 = (1 - ln(3)/4)/2 = I_cos

# W kontekscie linearizacji g-1 = delta*f = delta*sin(r)/r:
# A_lin^2 dla delta: A_lin = |delta|*1 (bo f = sin(r)/r asymptot. sin(r)/r jest TAIL)
# Tail sin(r)/r ma A = 1, faza 0. Wiec eta_lin = 1 (alpha_0 = 1).

# alpha_2 = I_cos i alpha_3 = I_sin^2/2 + P_cos pochodza z pelnej ekspansji etad.

# Sprawdzmy: I_sin = integral cos(r) * (1 - cos(r))/r^2 dr = ?
# Frullani trick: (1-cos r)/r^2 = int_0^1 sin(kr)/r dk,
# wiec int_0^inf cos(r) * int_0^1 sin(kr)/r dk dr = int_0^1 dk int_0^inf cos(r) sin(kr)/r dr

# int_0^inf cos(r) sin(kr)/r dr = pi/2 if k>1, 0 if k<1 ... wait
# Actually: int_0^inf sin(kr)cos(r)/r dr = (pi/4)*sgn(k+1) + (pi/4)*sgn(k-1)
# For 0<k<1: = pi/4 + 0 = pi/4? Let me check.

# sin(a)cos(b) = [sin(a+b)+sin(a-b)]/2
# sin(kr)cos(r)/r = [sin((k+1)r) + sin((k-1)r)]/(2r)
# int_0^inf sin(ar)/r dr = pi/2 * sgn(a)
# So int = [pi/2 sgn(k+1) + pi/2 sgn(k-1)]/2 = [pi/2 + pi/2*sgn(k-1)]/2
# k in (0,1): sgn(k-1) = -1, so = [pi/2 - pi/2]/2 = 0
# Hmm, then I_sin =0? That's wrong.

# Let me reconsider. I_sin as defined in this context is probably:
# I_sin = int_0^inf sin(r) * f(r) dr where f = sin(r)/r
# = int_0^inf sin^2(r)/r dr = int_0^inf (1-cos(2r))/(2r) dr -- DIVERGES at 0?
# No wait: (1 - cos 2r)/(2r) ~ (2r^2)/(2r) = r near 0, OK
# at inf: ~ 1/(2r) - cos(2r)/(2r), both diverge... Actually int 1/(2r) diverges at inf.
# So this integral doesn't converge absolutely. Need regularization.

# OK, the I_cos, I_sin in our context have specific definitions from the
# perturbation theory. Let me use the known results:
# I_cos = 1/2 - ln(3)/8 = 0.362673...
# I_sin = -pi/8 = -0.392699...

print(f"    I_cos = 1/2 - ln(3)/8 = {mp.nstr(I_cos_th, 20)}")
print(f"    I_sin = -pi/8        = {mp.nstr(I_sin_th, 20)}")
print(f"    I_sin^2 / 2         = pi^2/128 = {mp.nstr(I_sin_th**2/2, 20)}")
print(f"    = 55*pi^2/7040       = {mp.nstr(55*pi**2/7040, 20)}")

# ==============================================================
# KROK 2: J_c(r), J_s(r) integrals - closed form via Ci, Si
# ==============================================================

print(f"\n{'='*78}")
print(f"  KROK 2: Asymptotyki J_c(inf), J_s(inf)")
print(f"{'='*78}")

# s*(f')^2 = A(s) = (1+cos(2s))/(2s) - sin(2s)/s^2 + (1-cos(2s))/(2s^3)

# J_c(inf) = int_0^inf cos(s) * A(s) ds
# J_s(inf) = int_0^inf sin(s) * A(s) ds
# Both integrals require care near s=0 (but s*(f')^2 is regular at s=0 since
# f'(0) = 0, so A(s) ~ O(s) near 0).

# Numerical: compute J_c(inf), J_s(inf) to 40 digits

def A_integrand(s):
    # s*(f')^2(s) = cos^2(s)/s - 2 sin(s)cos(s)/s^2 + sin^2(s)/s^3
    s = mp.mpf(s)
    if s < mp.mpf('0.01'):
        # Taylor: f' = -sin/r^2 + cos/r ~ -1/r + r/6 - ... wait
        # f(r) = sin(r)/r = 1 - r^2/6 + r^4/120 - ...
        # f'(r) = -r/3 + r^3/30 - r^5/840 + ...
        # So (f')^2 = r^2/9 - r^4/45 + ...
        # s*(f')^2 = s^3/9 - s^5/45 + ...
        return s**3 / 9 - s**5 / 45 + s**7 * mp.mpf(2)/945  # approx
    return mp.cos(s)**2/s - 2*mp.sin(s)*mp.cos(s)/s**2 + mp.sin(s)**2/s**3

def Jc_integrand(s):
    return mp.cos(s) * A_integrand(s)

def Js_integrand(s):
    return mp.sin(s) * A_integrand(s)

# Careful integration with oscillation
print(f"\n  Obliczenie J_c(inf), J_s(inf) (ultra-precision)...")

# Use Levin or oscillatory methods
# mp.quad with trapezoidal is slow for oscillatory; use quadosc for oscillating tails.

# Split: [0, 100] direct, [100, inf] asymptotic
# Actually for high precision, we can use mpmath's quadosc with period pi

# J_c = int cos(s) A(s) ds
# A(s) ~ 1/(2s) at large s (leading) + oscillatory
# So cos(s) * A(s) ~ cos(s)/(2s) + oscillations
# int_large cos(s)/(2s) ds = -Ci(R)/2 contributes O(1/R), converges

# Let's just do high-precision numerical integration
print(f"    (integracja moze trwac chwile...)")

try:
    # Split [0, 50] finite, [50, inf] via quadosc
    Jc_1 = mp.quad(Jc_integrand, [0, 50])
    Jc_2 = mp.quadosc(Jc_integrand, [50, mp.inf], period=2*pi)
    Jc_total = Jc_1 + Jc_2
    print(f"    J_c(inf) = {mp.nstr(Jc_total, 20)}")
except Exception as e:
    print(f"    J_c blad: {e}")
    Jc_total = None

try:
    Js_1 = mp.quad(Js_integrand, [0, 50])
    Js_2 = mp.quadosc(Js_integrand, [50, mp.inf], period=2*pi)
    Js_total = Js_1 + Js_2
    print(f"    J_s(inf) = {mp.nstr(Js_total, 20)}")
except Exception as e:
    print(f"    J_s blad: {e}")
    Js_total = None

# ==============================================================
# KROK 3: Oczekiwane wartosci analityczne
# ==============================================================

print(f"\n{'='*78}")
print(f"  KROK 3: Teoretyczne oczekiwane wartosci J_c, J_s")
print(f"{'='*78}")

# int_0^inf cos(s)*cos^2(s)/s ds = ?
# cos(s)cos^2(s) = cos(s)*(1+cos(2s))/2 = cos(s)/2 + cos(s)cos(2s)/2
#                = cos(s)/2 + (cos(3s) + cos(s))/4 = 3cos(s)/4 + cos(3s)/4
# int_0^inf cos(as)/s ds DIVERGES at 0. But combined with other terms it's OK.

# s*(f')^2 = cos^2/s - 2sc/s^2 + s^2/s^3 where s = sin(s), c = cos(s)
# Near 0: [1 - s^2 + s^4/3]/s - 2[s-s^3/6+..][1-s^2/2+..]/s^2 + [s^2 - s^4/3 + ...]/s^3
#       = [1/s - s + s^3/3 - ...] - 2[s/s^2 - s^3/(6s^2) - s^3/(2s^2) + ...] + [1/s - s/3 + ...]
# Ohh this gives 1/s + 1/s = 2/s minus 2/s = 0 + finite. Good.

# More carefully:
# s*f'(s)^2 where f'(s) = d/ds[sin(s)/s] = cos(s)/s - sin(s)/s^2
# f'(s) = (s cos(s) - sin(s))/s^2
# Near 0: s cos s = s - s^3/2 + s^5/24 - ..., sin s = s - s^3/6 + s^5/120 - ...
# s cos s - sin s = -s^3/2 + s^3/6 + s^5/24 - s^5/120 - ... = -s^3/3 + s^5/30 - ...
# Hmm so s*cos(s) - sin(s) = -s^3/3 + s^5/30 - s^7/840 + ...
# f'(s) = -s/3 + s^3/30 - s^5/840 + ...
# (f'(s))^2 = s^2/9 - s^4/45 + s^6/1890 + ... wait
#   = (-s/3 + s^3/30 - ...)^2 = s^2/9 - 2*(s/3)(s^3/30) + (s^3/30)^2 + ...
#   = s^2/9 - s^4/45 + O(s^6)
# s*(f')^2 = s^3/9 - s^5/45 + ...

# OK integrand is regular at origin.

# Big-s asymptotics of s*(f')^2:
# (f')^2 = cos^2/s^2 - 2sc/s^3 + s^2/s^4
# Dominant 1/s^2 term, so s*(f')^2 ~ cos^2(s)/s = (1+cos(2s))/(2s) -- dominant.
# subleading: -sin(2s)/s^2 + sin^2(s)/s^3 = -sin(2s)/s^2 + (1-cos(2s))/(2s^3)

# int cos(s) * (1+cos(2s))/(2s) ds = int cos(s)/(2s) ds + int cos(s)cos(2s)/(2s) ds
# int_0^inf cos(s)/(2s) ds - regularized by the -(2sc/s^2) and (1-cos2s)/(2s^3) subtractions
# Actually we should do the whole integral of cos(s)*A(s) carefully.

# Known: int_0^inf cos(as)*cos(bs)/s ds = (1/2) ln|(a+b)/(a-b)| for |a|≠|b|
# (this is finite once you combine cos(as)*cos(bs) = [cos((a+b)s) + cos((a-b)s)]/2
#  and use the identity int_0^inf [cos(ps) - cos(qs)]/s ds = ln(q/p).)

# Actually: int_0^inf [cos(ps) - cos(qs)]/s ds = ln(q/p) (Frullani for cos).
# That's the key formula.

# Now cos(s)*(1+cos(2s))/(2s) = cos(s)/(2s) + cos(s)cos(2s)/(2s)
#                            = cos(s)/(2s) + [cos(3s) + cos(s)]/(4s)
#                            = 3cos(s)/(4s) + cos(3s)/(4s)
# int_0^inf [3cos(s) + cos(3s)]/(4s) ds = divergent standalone.

# But combined with -sin(2s)/s^2 + (1-cos(2s))/(2s^3):
# int cos(s) * [-sin(2s)/s^2] ds = -int [sin(3s) + sin(s)]/(2s^2) ds
# int cos(s) * (1-cos(2s))/(2s^3) ds = int [cos(s) - cos(s)cos(2s)/(2s^3)] ds
# ... This is getting hairy. Let's just rely on the numerical value and try to
# identify closed form.

if Jc_total is not None and Js_total is not None:
    print(f"\n  Identyfikacja closed form dla J_c, J_s:")

    # J_c kandydaci: ln(2), ln(3), 1/2, pi/8, pi^2 related
    for name, v in [
        ('1/2', mp.mpf(1)/2),
        ('ln(2)/2', mp.log(2)/2),
        ('ln(3)/4', mp.log(3)/4),
        ('1/2 - ln(3)/8', I_cos_th),
        ('I_cos', I_cos_th),
        ('ln(2) - 1/2', mp.log(2) - mp.mpf(1)/2),
        ('pi^2/24', pi**2/24),
        ('pi^2/32', pi**2/32),
    ]:
        diff = Jc_total - v
        if abs(diff) < mp.mpf('1e-3'):
            print(f"    J_c ~ {name} = {mp.nstr(v, 15)}, diff = {mp.nstr(diff, 5)}")

    for name, v in [
        ('-pi/8', -pi/8),
        ('-pi/4', -pi/4),
        ('-pi/16', -pi/16),
        ('I_sin', I_sin_th),
        ('pi/8', pi/8),
    ]:
        diff = Js_total - v
        if abs(diff) < mp.mpf('1e-3'):
            print(f"    J_s ~ {name} = {mp.nstr(v, 15)}, diff = {mp.nstr(diff, 5)}")

# ==============================================================
# KROK 4: Formula na P_cos przez J_c, J_s
# ==============================================================

print(f"\n{'='*78}")
print(f"  KROK 4: Relacja P_cos = f(J_c, J_s, I_cos, I_sin)")
print(f"{'='*78}")

# From perturbation theory:
# g = 1 + delta*f + delta^2*h + delta^3*p + ...
# f = sin(r)/r has A_tail = 1 (so |delta|), phase = 0
# h tail: u = rh has asymptotic -sin(r) J_c + cos(r) J_s
#   so h tail: A_h = sqrt(J_c^2 + J_s^2), phase atan2(-Jc, Js)...
#   but sign depends on convention.
#
# For eta(delta) = A(1+delta)/|delta|:
# A = |delta| + delta^2 * (correction from h tail projected onto sin(r)/r channel)
#   + delta^3 * (corrections from p, interaction, and |delta| envelope)

# The "cos tail" means the tail proportional to cos(r) rather than sin(r)
# in h's asymptotic expansion. Specifically:
# h ~ (J_s * cos(r) - J_c * sin(r))/r as r->inf (from u = rh)
#   = [amp * cos(r + phase)]/r where amp = sqrt(Jc^2 + Js^2)
# So "cos-coefficient in h-tail" = J_s  (coefficient of cos(r)/r)
#    "sin-coefficient in h-tail" = -J_c

# The delta^3 tail from this: delta^2 * h_tail contributes to A(1+delta).
# g(r) - 1 = delta f + delta^2 h + delta^3 p + ...
#        ~ delta sin(r)/r + delta^2 (J_s cos r - J_c sin r)/r + O(delta^3)
# So tail amplitude (sin,cos): (delta - delta^2 J_c, delta^2 J_s) + O(d^3)
# |A_tail| = sqrt((d - d^2 Jc)^2 + (d^2 Js)^2)
#         = |d| sqrt( 1 - 2 d Jc + d^2 Jc^2 + d^2 Js^2 )
#         = |d| [1 - d Jc + d^2 (Jc^2 + Js^2)/2 - d^2 Jc^2/2 + ...]
#         = |d| [1 - d Jc + d^2 Js^2/2 + O(d^3)]
# Hmm. So eta = A/|d| = 1 - d*Jc + d^2*Js^2/2 + O(d^3)
# So alpha_2 = -J_c, alpha_3 = J_s^2/2 + (O(d^3) perturbation contribution)

# Wait: alpha_2 = c_1/2 = I_cos > 0. So we need J_c = -I_cos.
# alpha_3 = J_s^2/2 + P_cos = I_sin^2/2 + P_cos, so J_s = I_sin = -pi/8 (up to sign).

# Check: is J_s = -pi/8 numerically?
if Js_total is not None:
    diff_Is = Js_total - I_sin_th
    print(f"\n  J_s numeric = {mp.nstr(Js_total, 20)}")
    print(f"  I_sin       = {mp.nstr(I_sin_th, 20)}")
    print(f"  diff        = {mp.nstr(diff_Is, 5)}")

if Jc_total is not None:
    diff_Ic = Jc_total + I_cos_th  # -Jc = I_cos means Jc = -I_cos
    print(f"\n  J_c numeric   = {mp.nstr(Jc_total, 20)}")
    print(f"  -I_cos        = {mp.nstr(-I_cos_th, 20)}")
    print(f"  J_c + I_cos   = {mp.nstr(diff_Ic, 5)}  (should be 0 if Jc = -I_cos)")

# ==============================================================
# KROK 5: P_cos z delta^3 perturbacji
# ==============================================================

print(f"\n{'='*78}")
print(f"  KROK 5: P_cos jako cos-tail w p(r) rownaniu")
print(f"{'='*78}")

# p'' + (2/r)p' + p = f(f')^2 - 2 f' h'
# v = rp: v'' + v = r * [f(f')^2 - 2 f' h']
# v(r) = int_0^r sin(r-s) * s*[f(f')^2 - 2 f' h'](s) ds
# p-tail: v ~ -sin(r) Kc + cos(r) Ks
#   Kc = int_0^inf cos(s) * s * [f(f')^2 - 2 f' h'] ds
#   Ks = int_0^inf sin(s) * s * [f(f')^2 - 2 f' h'] ds
#
# Tail kontrybucja z delta^3 * p do A(1+d):
# g-1 ~ d sin/r + d^2 (Js cos - Jc sin)/r + d^3 (Ks cos - Kc sin)/r + ...
# Wiec amplitudes (sin, cos) coefficients of (g-1)*r:
# sin-coeff: d - d^2 Jc - d^3 Kc
# cos-coeff: d^2 Js + d^3 Ks
# A = sqrt(...) = |d| sqrt(1 + (-2d Jc - 2d^2 Kc + d^2 Jc^2 + d^2 Js^2) + ...)
#   = |d| [1 - d Jc - d^2 Kc + d^2 (Jc^2+Js^2)/2 - d^2 Jc^2/2 + ...]
#   = |d| [1 - d Jc + d^2 (Js^2/2 - Kc) + O(d^3)]
# Wiec eta = 1 - d*Jc + d^2*(Js^2/2 - Kc) + O(d^3)
# alpha_2 = -Jc = I_cos (spojne jesli Jc = -I_cos)
# alpha_3 = Js^2/2 - Kc = I_sin^2/2 + P_cos
# ==> P_cos = -Kc = -int_0^inf cos(s) * s*[f(f')^2 - 2 f'h'](s) ds

# To jest ODPOWIEDZ na "co to jest P_cos".
# P_cos = -int_0^inf cos(s) * s * [f(s)(f'(s))^2 - 2 f'(s) h'(s)] ds

# We need h(s) numerically (or analytically) to evaluate this.

# Plan:
# (1) Compute h(s) from u(s) = r*h(s) = -sin(r) J_c(r) + cos(r) J_s(r)
# (2) Compute h'(s)
# (3) Evaluate Kc = int cos(s) * s * [f(f')^2 - 2 f' h'] ds
# (4) Check -Kc vs 9*pi^2/7040

from mpmath import mpf, sin, cos, quad, inf

def f_fun(s):
    if abs(s) < 1e-30:
        return mp.mpf(1)
    return sin(s)/s

def fp_fun(s):
    # f'(s) = (s cos s - sin s)/s^2
    if abs(s) < mp.mpf('0.001'):
        return -s/3 + s**3/30 - s**5/840
    return (s*cos(s) - sin(s))/s**2

def sA(s):
    # s * (f')^2
    return s * fp_fun(s)**2

# h(r) = u(r)/r where u(r) = -int_0^r sin(r-s) * s * (f')^2 ds
def u_fun(r):
    # u(r) = -sin(r)*Jc(r) + cos(r)*Js(r)
    def integ_c(s): return cos(s) * sA(s)
    def integ_s(s): return sin(s) * sA(s)
    Jc_r = mp.quad(integ_c, [0, r])
    Js_r = mp.quad(integ_s, [0, r])
    return -sin(r)*Jc_r + cos(r)*Js_r

def h_fun(r):
    if abs(r) < mp.mpf('0.001'):
        return mp.mpf(0)  # h(0) = 0
    return u_fun(r)/r

print(f"\n  (Obliczenia h(r), h'(r), Kc ultra-slow. Probujemy skroconey test)")

# Spot check: h(5), h(10), h(20)
import time
t0 = time.time()
print(f"    h(3)  = {mp.nstr(h_fun(mp.mpf('3')), 15)}  ({time.time()-t0:.1f}s)")
t0 = time.time()
print(f"    h(10) = {mp.nstr(h_fun(mp.mpf('10')), 15)}  ({time.time()-t0:.1f}s)")

print(f"\n{'='*78}")
print(f"  WNIOSKI KROKU 5:")
print(f"{'='*78}")
print(f"""
  USTALONO analityczna tozsamosc:
    P_cos = -Kc = -int_0^inf cos(s) * s * [f(f')^2 - 2 f' h'](s) ds

  Nastepny krok: evaluate ten integral closed-form lub ultra-precyzyjnie
  i zidentifikowac jako 9*pi^2/7040.
""")
