#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c2_refined_precision.py
===========================
Wysokoprecyzyjny pomiar alpha_3 z kontrola bias.

STRATEGIA (oparta na obserwacji z r6_c2_identification_scan.py):
  sym/delta^2 plateau-je na ~0.08971 przy delta in [0.003, 0.012],
  a Richardson rozjezdza sie do 0.08994 bo uzyl noisy-ch malych delta (<0.003).

POPRAWIONA EKSTRAKCJA alpha_3:
  1. Pomiar na "dobrym" oknie delta in [0.003, 0.012] (malo szumu, male O(delta^4))
  2. Zmienny r_max (400, 600, 800) -> sprawdzic bias finite-r_max
  3. Rozszerzony tail fit: u(r) = c1*cos(r) + c2*sin(r) + d1*cos(r)/r + d2*sin(r)/r
     (uwzglednia O(1/r^2) korekcje do A_tail)
  4. Porownanie: czy alpha_3 = pi^2/110 = 0.0897237, czy inny?

Docelowa precyzja: bias < 3e-5 (rozroznienie pi^2/110 od sasiednich kandydatow)

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
    """Standardowy fit u(r)=(g-1)*r = c1*cos + c2*sin"""
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max_win)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    return math.sqrt(coef[0]**2 + coef[1]**2)

def extract_atail_extended(sol, r_min, r_max_win):
    """Rozszerzony fit uwzgledniajacy O(1/r) korekcje:
       u(r) = c1*cos(r) + c2*sin(r) + (d1*cos(r) + d2*sin(r))/r
       Amplituda A = sqrt(c1^2 + c2^2)"""
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max_win)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([
        np.cos(r_f), np.sin(r_f),
        np.cos(r_f)/r_f, np.sin(r_f)/r_f
    ])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    return math.sqrt(coef[0]**2 + coef[1]**2)

def measure_eta(delta, r_max, n_pts, extractor):
    sol = solve_ode(1.0 + delta, r_max, n_pts)
    if not sol.success: return None
    A = extractor(sol, 100.0, r_max - 50.0)
    return A / abs(delta)

print("=" * 80)
print("  alpha_3: WYSOKA PRECYZJA z kontrola bias (multiple r_max)")
print("=" * 80)

# Wybieramy "dobre" delta: nie za male (zeby nie miec szumu w y = sym-1 ~ 1e-7)
# ale tez nie za duze (zeby O(delta^4) byl maly)
deltas = np.array([0.003, 0.004, 0.005, 0.006, 0.008, 0.010, 0.012])

c1_theory = 1.0 - math.log(3.0)/4.0

# Rozne r_max -> Richardson na r_max -> 0
scenarios = [
    ("r_max=400,  basic", 400.0,  160000, extract_atail_basic),
    ("r_max=600,  basic", 600.0,  240000, extract_atail_basic),
    ("r_max=800,  basic", 800.0,  320000, extract_atail_basic),
    ("r_max=400,  ext  ", 400.0,  160000, extract_atail_extended),
    ("r_max=600,  ext  ", 600.0,  240000, extract_atail_extended),
    ("r_max=800,  ext  ", 800.0,  320000, extract_atail_extended),
]

results = {}
for label, r_max, n_pts, extractor in scenarios:
    print(f"\n[{label}]")
    sym_minus_1 = []
    eta_p_list = []
    eta_m_list = []
    for d in deltas:
        e_p = measure_eta(+d, r_max, n_pts, extractor)
        e_m = measure_eta(-d, r_max, n_pts, extractor)
        if e_p is None or e_m is None:
            continue
        sym = (e_p + e_m) / 2.0
        sym_minus_1.append(sym - 1.0)
        eta_p_list.append(e_p)
        eta_m_list.append(e_m)

    y = np.array(sym_minus_1)
    d_arr = deltas
    # fit y = a3*d^2 + a5*d^4
    X = np.column_stack([d_arr**2, d_arr**4])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a3 = coef[0]
    a5 = coef[1]
    results[label] = (a3, a5, eta_p_list, eta_m_list)

    # alpha_2 verification: asym slope
    asym_slopes = [(eta_p_list[i] - eta_m_list[i])/(2*d_arr[i]) for i in range(len(d_arr))]
    alpha_2 = np.mean(asym_slopes)
    print(f"  alpha_2 (asym slope mean) = {alpha_2:.10f}  (theory c1/2 = {c1_theory/2:.10f})")
    print(f"  diff alpha_2 - c1/2      = {alpha_2 - c1_theory/2:+.3e}")
    print(f"  alpha_3 (poly d^2 + d^4)  = {a3:.10f}")
    print(f"  alpha_5                   = {a5:+.6f}")

print(f"\n{'='*80}")
print(f"  KANDYDACI (z kontrola bias):")
print(f"{'='*80}")
alpha3_vals = [results[s[0]][0] for s in scenarios]
print(f"\n  {'Scenariusz':25s}  {'alpha_3':>14s}")
for s in scenarios:
    a3 = results[s[0]][0]
    print(f"  {s[0]:25s}  {a3:14.10f}")

# Richardson in r_max
# Expected: alpha_3(r_max) = alpha_3_infinity + k/r_max^2 (or other power)
# Compare basic at 400, 600, 800:
a3_400 = results["r_max=400,  basic"][0]
a3_600 = results["r_max=600,  basic"][0]
a3_800 = results["r_max=800,  basic"][0]
print(f"\n  Richardson r_max -> inf (basic):")
print(f"    r_max=400: {a3_400:.10f}")
print(f"    r_max=600: {a3_600:.10f}")
print(f"    r_max=800: {a3_800:.10f}")

# Fit alpha_3(r_max) = A + B / r_max^p, try p in [1, 2, 3]
for p in [1, 2, 3]:
    rm_arr = np.array([400, 600, 800])
    a3_arr = np.array([a3_400, a3_600, a3_800])
    X = np.column_stack([np.ones(3), 1.0/rm_arr**p])
    coef, *_ = np.linalg.lstsq(X, a3_arr, rcond=None)
    print(f"    fit alpha_3 = A + B/r_max^{p}: A = {coef[0]:.10f}, B = {coef[1]:+.6f}")

# Same for extended
a3_400_e = results["r_max=400,  ext  "][0]
a3_600_e = results["r_max=600,  ext  "][0]
a3_800_e = results["r_max=800,  ext  "][0]
print(f"\n  Richardson r_max -> inf (extended fit):")
print(f"    r_max=400: {a3_400_e:.10f}")
print(f"    r_max=600: {a3_600_e:.10f}")
print(f"    r_max=800: {a3_800_e:.10f}")

for p in [1, 2, 3]:
    rm_arr = np.array([400, 600, 800])
    a3_arr = np.array([a3_400_e, a3_600_e, a3_800_e])
    X = np.column_stack([np.ones(3), 1.0/rm_arr**p])
    coef, *_ = np.linalg.lstsq(X, a3_arr, rcond=None)
    print(f"    fit alpha_3 = A + B/r_max^{p}: A = {coef[0]:.10f}, B = {coef[1]:+.6f}")

# Porownanie z kandydatami
print(f"\n{'='*80}")
print(f"  KANDYDACI CLOSED-FORM (porownanie z najlepszym pomiarem):")
print(f"{'='*80}")
ln3 = math.log(3)
cands = {
    'pi^2/110':              math.pi**2 / 110,
    '(1-c1)/12 + 8(1-c1)^2/9': (1-c1_theory)/12 + 8*(1-c1_theory)**2/9,
    'ln(3)/48 + ln(3)^2/18':   ln3/48 + ln3**2/18,
    '(pi-1)/24':             (math.pi-1)/24,
    '(2-ln3)/10':            (2-ln3)/10,
    'pi/35':                 math.pi / 35,
    'c1^2/6':                c1_theory**2 / 6,
}
best = a3_800_e  # assume r_max=800 with extended fit is most accurate
print(f"\n  Best estimate (r_max=800, ext): {best:.10f}")
print(f"\n  {'Kandydat':35s}  {'Wartosc':>14s}  {'Diff':>12s}  {'Rel':>10s}")
print(f"  {'-'*35}  {'-'*14}  {'-'*12}  {'-'*10}")
for name, v in sorted(cands.items(), key=lambda kv: abs(kv[1]-best)):
    d = v - best
    rel = abs(d) / abs(best)
    print(f"  {name:35s}  {v:14.10f}  {d:+12.3e}  {rel:10.2e}")
