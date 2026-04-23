#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c2_ultra_precision.py
=========================
Ultra-precyzyjny pomiar alpha_3 do rozstrzygniecia: czy alpha_3 = pi^2/110?

Z r6_c2_refined_precision.py:
  alpha_3 (basic, r_max=400): 0.0897100596
  alpha_3 (basic, r_max=600): 0.0897146689
  alpha_3 (basic, r_max=800): 0.0897170601
  Richardson A + B/r_max^2:   A = 0.0897190722
  pi^2/110                  = 0.0897236764   (diff 5e-6 od A)

Dodajemy r_max=1200, 1600 -> powinno bardziej precyzyjnie oszacowac A_inf.

Cel: rozstrzygniecie czy alpha_3 = pi^2/110 w granicach <1e-5
LUB zobaczyc czy A_inf odbiega od pi^2/110.

Author: Claudian
Date: 2026-04-16
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp

def rhs(r, y):
    g, gp = y
    if g < 1e-12: g = 1e-12
    if r < 1e-13:
        gpp = (1.0 - g) / 4.0
    else:
        gpp = (1.0 - g) - (1.0/g) * gp**2 - (2.0/r) * gp
    return [gp, gpp]

def solve_ode(g0, r_max, n_pts):
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=1e-13, atol=1e-15, max_step=0.015)
    return sol

def extract_atail_basic(sol, r_min, r_max_win):
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max_win)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    return math.sqrt(coef[0]**2 + coef[1]**2)

def measure_eta(delta, r_max, n_pts):
    sol = solve_ode(1.0 + delta, r_max, n_pts)
    if not sol.success: return None
    A = extract_atail_basic(sol, 100.0, r_max - 50.0)
    return A / abs(delta)

print("=" * 78)
print("  alpha_3: ULTRA PRECYZJA (r_max do 1600)")
print("=" * 78)

deltas = np.array([0.003, 0.004, 0.005, 0.006, 0.008, 0.010, 0.012])
c1_th = 1.0 - math.log(3.0)/4.0

r_max_list = [400, 600, 800, 1200, 1600]
# n_pts scaled so r-spacing is ~constant
n_pts_list = [int(400*rm) for rm in r_max_list]  # ~= constant 400 pts/unit

alpha3_dict = {}
alpha2_dict = {}
for r_max, n_pts in zip(r_max_list, n_pts_list):
    print(f"\n  [r_max={r_max}, n_pts={n_pts}]")
    sym_arr = []
    eta_p_arr = []
    eta_m_arr = []
    for d in deltas:
        e_p = measure_eta(+d, r_max, n_pts)
        e_m = measure_eta(-d, r_max, n_pts)
        sym = (e_p + e_m) / 2.0
        sym_arr.append(sym - 1.0)
        eta_p_arr.append(e_p)
        eta_m_arr.append(e_m)
    y = np.array(sym_arr)
    X = np.column_stack([deltas**2, deltas**4])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a3 = coef[0]
    asym_slopes = [(eta_p_arr[i] - eta_m_arr[i])/(2*deltas[i]) for i in range(len(deltas))]
    alpha_2 = np.mean(asym_slopes)
    alpha3_dict[r_max] = a3
    alpha2_dict[r_max] = alpha_2
    print(f"    alpha_2 = {alpha_2:.12f}  (diff c1/2 = {alpha_2 - c1_th/2:+.3e})")
    print(f"    alpha_3 = {a3:.12f}")

print(f"\n{'='*78}")
print(f"  RICHARDSON EXTRAPOLATION r_max -> inf")
print(f"{'='*78}")

rm_arr = np.array(r_max_list, dtype=float)
a3_arr = np.array([alpha3_dict[r] for r in r_max_list])
a2_arr = np.array([alpha2_dict[r] for r in r_max_list])

print(f"\n  Zebrane wartosci:")
print(f"  {'r_max':>8s}  {'alpha_2':>14s}  {'alpha_3':>14s}  {'diff a2':>12s}  {'diff a3(pi^2/110)':>18s}")
for rm in r_max_list:
    a2 = alpha2_dict[rm]
    a3 = alpha3_dict[rm]
    print(f"  {rm:8.0f}  {a2:14.10f}  {a3:14.10f}  {a2-c1_th/2:+12.3e}  {a3 - math.pi**2/110:+18.3e}")

# Fit alpha_3 = A + B/r_max^p
for p in [1, 2, 3, 4]:
    X = np.column_stack([np.ones(len(rm_arr)), 1.0/rm_arr**p])
    coef, res, *_ = np.linalg.lstsq(X, a3_arr, rcond=None)
    # residual RMS
    pred = X @ coef
    rms = np.sqrt(np.mean((a3_arr - pred)**2))
    print(f"\n  fit alpha_3(r_max) = A + B/r_max^{p}: A = {coef[0]:.12f}, B = {coef[1]:+.6e}, RMS = {rms:.3e}")

# Same dla alpha_2
print(f"\n  === alpha_2 extrapolation (weryfikacja, powinno = c1/2 = {c1_th/2:.12f}) ===")
for p in [1, 2, 3]:
    X = np.column_stack([np.ones(len(rm_arr)), 1.0/rm_arr**p])
    coef, *_ = np.linalg.lstsq(X, a2_arr, rcond=None)
    print(f"  fit alpha_2 = A + B/r_max^{p}: A = {coef[0]:.12f}, diff c1/2 = {coef[0]-c1_th/2:+.3e}")

# Compare final A with candidates
print(f"\n{'='*78}")
print(f"  FINAL CANDIDATE COMPARISON")
print(f"{'='*78}")

# Best alpha_3 estimate: use 1/r_max^2 fit (motivated by O(1/r^2) finite-r correction)
X = np.column_stack([np.ones(len(rm_arr)), 1.0/rm_arr**2])
coef, *_ = np.linalg.lstsq(X, a3_arr, rcond=None)
alpha3_best = coef[0]

ln3 = math.log(3)
cands = {
    'pi^2/110':                math.pi**2 / 110,
    'pi/35':                   math.pi / 35,
    'ln(3)/48 + ln(3)^2/18':   ln3/48 + ln3**2/18,
    '(1-c1)/12 + 8(1-c1)^2/9': (1-c1_th)/12 + 8*(1-c1_th)**2/9,
    '(pi-1)/24':               (math.pi-1)/24,
    '1/(4*pi) + 1/100':        1/(4*math.pi) + 1/100,
    '(3-ln3)^2/43':            (3-ln3)**2 / 43,
    'c1^2 * (9/20) * (1-c1)':  c1_th**2 * (9/20) * (1-c1_th),
    '(9/20)*c1*(1-c1)':        (9.0/20)*c1_th*(1-c1_th),
    '(3 - e)/3':               (3-math.e)/3,
    'pi^2/(11*10)':            math.pi**2/(11*10),  # same as pi^2/110
    '1/(4*pi) + 1/100':        1/(4*math.pi) + 0.01,
    '(pi/2)^2/27.5':           (math.pi/2)**2/27.5,
    '(9/(pi+7))^2':            (9/(math.pi+7))**2,
}

print(f"\n  Best alpha_3 estimate (extrap 1/r^2): {alpha3_best:.12f}")
print(f"\n  {'Kandydat':35s}  {'Wartosc':>14s}  {'Diff':>12s}  {'Rel':>10s}")
print(f"  {'-'*35}  {'-'*14}  {'-'*12}  {'-'*10}")
for name, v in sorted(set(cands.items()), key=lambda kv: abs(kv[1]-alpha3_best)):
    d = v - alpha3_best
    rel = abs(d) / abs(alpha3_best)
    marker = " <--" if abs(d) < 5e-5 else ""
    print(f"  {name:35s}  {v:14.10f}  {d:+12.3e}  {rel:10.2e}{marker}")

# Precyzja a2 powinna validate cala procedure
X = np.column_stack([np.ones(len(rm_arr)), 1.0/rm_arr**2])
coef2, *_ = np.linalg.lstsq(X, a2_arr, rcond=None)
alpha2_best = coef2[0]
print(f"\n  Alpha_2 best = {alpha2_best:.12f}")
print(f"  c1/2         = {c1_th/2:.12f}")
print(f"  Diff         = {alpha2_best - c1_th/2:+.3e}")
print(f"\n  Jesli alpha_2 zgadza sie z c1/2 do < 1e-5, to metoda extrapolacji jest")
print(f"  godna zaufania i alpha_3 wynik tez bedzie do < 1e-5 rzetelny.")
