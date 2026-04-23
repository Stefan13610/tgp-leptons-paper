#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c3_perturbative_h_and_p.py
==============================

Analityczne wyprowadzenie α_3 z perturbacji ODE.

PLAN:
  1. Rozwiazac h(r) analitycznie (lub numerycznie z wysoka precyzja)
     h'' + (2/r)h' + h = -(f')^2,  f = sin(r)/r
  2. Zweryfikowac asymptoty: u_h(r) = r*h(r) ~ sin(r)*I_cos - cos(r)*I_sin
     gdzie I_cos = 1/2 - ln(3)/8 i I_sin = -pi/8
  3. Rozwiazac p(r) z O(delta^3):
     p'' + (2/r)p' + p = -2*f'*h' + f*(f')^2
     u_p(r) = r*p(r) ~ sin(r)*P_cos - cos(r)*P_sin
  4. Porownac I_sin^2/2 + P_cos z alpha_3 = pi^2/110 = 0.0897237

Jesli zgodnosc do 10^-6, to:
  - Potwierdzamy strukture perturbacyjna
  - Mozemy szukac ZAMKNIETEJ formy dla P_cos (9*pi^2/7040?)

Author: Claudian
Date: 2026-04-16
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp, quad
from scipy.special import sici  # Si(x), Ci(x) integrals

# -----------------------------------------------------------------
# ANALITYCZNE KONTROLE
# -----------------------------------------------------------------
c1_th = 1.0 - math.log(3.0)/4.0  # 0.72534693
I_cos_th = 0.5 - math.log(3.0)/8.0  # 0.36267346 (proven in c1 breakthrough)
I_sin_th = -math.pi / 8.0           # -0.39269908 (proven in c1 breakthrough)
alpha3_target = math.pi**2 / 110.0  # 0.08972368

print("=" * 78)
print("  DERYWACJA alpha_3 z perturbacji O(delta^2) + O(delta^3)")
print("=" * 78)
print(f"\n  STALE ANALITYCZNE (z breakthrough c1):")
print(f"    I_cos = 1/2 - ln(3)/8 = {I_cos_th:.15f}")
print(f"    I_sin = -pi/8         = {I_sin_th:.15f}")
print(f"    c1 = 2*I_cos          = {2*I_cos_th:.15f}")
print(f"    I_sin^2/2             = {I_sin_th**2/2:.15f}  (pi^2/128)")
print(f"    alpha_3 target (pi^2/110) = {alpha3_target:.15f}")
print(f"    P_cos target = alpha_3 - I_sin^2/2 = {alpha3_target - I_sin_th**2/2:.15f}")
print(f"    (= 9*pi^2/7040 = {9*math.pi**2/7040:.15f} - sprawdzic)")

# -----------------------------------------------------------------
# 1. SOLVE h(r) NUMERYCZNIE
# -----------------------------------------------------------------
# f(r) = sin(r)/r,  f'(r) = cos(r)/r - sin(r)/r^2
# h'' + (2/r) h' + h = -(f')^2
# Substitucja u = r*h -> u'' + u = -r*(f')^2
# u(0) = 0, u'(0) = 0  (h regularne przy r=0, h(0) = skonczona)
# Ale trzeba sprawdzic: f'(0) = lim (cos(r)/r - sin(r)/r^2) = ?
# Rozwiniecie f(r) = 1 - r^2/6 + r^4/120 - ...
# f'(r) = -r/3 + r^3/30 - ... -> f'(0) = 0
# (f')^2 at r=0: 0
# So source near r=0: -r*(f')^2 ~ -r*(r/3)^2 = -r^3/9 (regular)

def source_u(r):
    """Zrodlo dla u'' + u = q(r), q(r) = -r*(f'(r))^2"""
    if r < 1e-8:
        return -r**3 / 9.0  # rozwiniecie Taylora
    fp = math.cos(r)/r - math.sin(r)/r**2
    return -r * fp * fp

def rhs_u(r, y):
    u, up = y
    return [up, -u + source_u(r)]

r_max = 2000.0
n_pts = 400000
r_eval = np.linspace(1e-8, r_max, n_pts)
print(f"\n  Rozwiazuje u'' + u = -r*(f')^2, r_max={r_max}")
sol_u = solve_ivp(rhs_u, (1e-8, r_max), [0.0, 0.0],
                  method='DOP853', t_eval=r_eval,
                  rtol=1e-12, atol=1e-15, max_step=0.01)
u = sol_u.y[0]
up = sol_u.y[1]

# Asymptotyka u(r) -> sin(r)*I_cos - cos(r)*I_sin
# Fit u = A*sin + B*cos w oknie [500, 1900]
mask = (sol_u.t >= 500.0) & (sol_u.t <= 1900.0)
r_f = sol_u.t[mask]
u_f = u[mask]
X = np.column_stack([np.sin(r_f), np.cos(r_f)])
coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
I_cos_num = coef[0]
I_sin_num = -coef[1]  # bo X drugi kolumna to cos, a fit to -cos * I_sin
print(f"\n  Weryfikacja: u(r) = sin(r)*I_cos - cos(r)*I_sin")
print(f"    I_cos (numer.) = {I_cos_num:.12f}")
print(f"    I_cos (theor.) = {I_cos_th:.12f}  diff = {I_cos_num - I_cos_th:+.3e}")
print(f"    I_sin (numer.) = {I_sin_num:.12f}")
print(f"    I_sin (theor.) = {I_sin_th:.12f}  diff = {I_sin_num - I_sin_th:+.3e}")

# Ekstrakcja h(r) = u(r)/r, h'(r) z up:
# h = u/r,  h' = (u' - h)/r = (up - u/r)/r = (up*r - u)/r^2
def h_of_r(r, u, up):
    return u / r
def hp_of_r(r, u, up):
    return (up * r - u) / r**2

# -----------------------------------------------------------------
# 2. SOLVE p(r) NUMERYCZNIE
# -----------------------------------------------------------------
# p'' + (2/r)p' + p = -2*f'*h' + f*(f')^2
# v = r*p -> v'' + v = r*[-2*f'*h' + f*(f')^2]
# v(0) = 0, v'(0) = 0
# Interpolacja u, up do obliczenia h, hp przy dowolnym r

from scipy.interpolate import CubicSpline

spline_u = CubicSpline(sol_u.t, u)
spline_up = CubicSpline(sol_u.t, up)

def source_v(r):
    """Zrodlo dla v'' + v = r*Q(r),  Q = -2*f'*h' + f*(f')^2"""
    if r < 1e-7:
        return 0.0  # regularne
    f_val = math.sin(r)/r
    fp = math.cos(r)/r - math.sin(r)/r**2
    u_val = float(spline_u(r))
    up_val = float(spline_up(r))
    h_val = u_val / r
    hp_val = (up_val * r - u_val) / r**2
    Q = -2.0 * fp * hp_val + f_val * fp * fp
    return r * Q

def rhs_v(r, y):
    v, vp = y
    return [vp, -v + source_v(r)]

# Rozwiazuje tylko w oknie [eps, r_max] gdzie mamy spline dla u
r_max_v = 1800.0  # zostaw margines dla splinu
r_eval_v = np.linspace(1e-7, r_max_v, 200000)
print(f"\n  Rozwiazuje v'' + v = r*(-2*f'*h' + f*(f')^2), r_max={r_max_v}")
sol_v = solve_ivp(rhs_v, (1e-7, r_max_v), [0.0, 0.0],
                  method='DOP853', t_eval=r_eval_v,
                  rtol=1e-11, atol=1e-14, max_step=0.02)
v = sol_v.y[0]
vp = sol_v.y[1]

# Asymptotyka v(r) -> sin(r)*P_cos - cos(r)*P_sin
mask_v = (sol_v.t >= 500.0) & (sol_v.t <= 1700.0)
r_fv = sol_v.t[mask_v]
v_f = v[mask_v]
Xv = np.column_stack([np.sin(r_fv), np.cos(r_fv)])
coef_v, *_ = np.linalg.lstsq(Xv, v_f, rcond=None)
P_cos_num = coef_v[0]
P_sin_num = -coef_v[1]
print(f"\n  Asymptotyka v(r) = sin(r)*P_cos - cos(r)*P_sin:")
print(f"    P_cos (numer.) = {P_cos_num:.12f}")
print(f"    P_sin (numer.) = {P_sin_num:.12f}")

# -----------------------------------------------------------------
# 3. WERYFIKACJA alpha_3 = I_sin^2/2 + P_cos
# -----------------------------------------------------------------
alpha3_predicted = I_sin_th**2 / 2 + P_cos_num
print(f"\n{'='*78}")
print(f"  WERYFIKACJA alpha_3 = I_sin^2/2 + P_cos")
print(f"{'='*78}")
print(f"\n  I_sin^2/2       = {I_sin_th**2/2:.12f}")
print(f"  P_cos (numer.)  = {P_cos_num:.12f}")
print(f"  SUM             = {alpha3_predicted:.12f}")
print(f"  alpha_3 target  = {alpha3_target:.12f}  (pi^2/110)")
print(f"  diff            = {alpha3_predicted - alpha3_target:+.6e}")

# Tez sprawdz przez hipoteze P_cos = 9*pi^2/7040
P_cos_hypo = 9 * math.pi**2 / 7040
print(f"\n  Hipoteza 9*pi^2/7040:")
print(f"    P_cos (numer.)  = {P_cos_num:.12f}")
print(f"    9*pi^2/7040     = {P_cos_hypo:.12f}")
print(f"    diff            = {P_cos_num - P_cos_hypo:+.3e}")

# Inne kandydaty dla P_cos
print(f"\n  Inne kandydaty P_cos:")
cands = {
    '9*pi^2/7040':         9*math.pi**2/7040,
    'pi^2/(8*99)':         math.pi**2/(8*99),
    'pi^2/782':            math.pi**2/782,
    '1/(8*pi^2) ~ 0.01267': 1.0/(8*math.pi**2),
    '(pi/2)^2/196':        (math.pi/2)**2/196,
    'pi^2/110 - pi^2/128': math.pi**2*(1/110 - 1/128),
    'ln(2)/55':            math.log(2)/55,
    'ln(3)/87':            math.log(3)/87,
    'ln(3)^2 / (64+32)':   math.log(3)**2/96,
    '1/(6*c1^2*16)':       1/(6*c1_th**2*16),
}
for name, v in sorted(cands.items(), key=lambda kv: abs(kv[1]-P_cos_num)):
    diff = v - P_cos_num
    print(f"    {name:30s} = {v:.10f}  diff = {diff:+.3e}")
