#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c2_identification_scan.py
==============================

Identyfikacja closed-form dla alpha_3 = 0.089710 (pomiar z r6_c2_perturbative_coeffs.py).

Strategia:
  1. Zmierz alpha_3 z jeszcze wyzsza precyzja (mniejsze delta, wiecej punktow)
  2. Dokladniejszy Richardson extrapolation delta -> 0
  3. Skan kandydatow kombinacji: rationals * {pi^k, ln(2), ln(3), Catalan, zeta(3)}
  4. Zidentyfikowac czy alpha_3 ma regulacje typu c1^n * rational + pi/ln term

Wsparcie na PSLQ-styl scan (nie mamy mpmath PSLQ ale rozny scan kombinacji).

Author: Claudian
Date: 2026-04-16
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp
from itertools import product

# ==================== Re-use solver ====================
def rhs(r, y):
    g, gp = y
    if g < 1e-12: g = 1e-12
    if r < 1e-13:
        gpp = (1.0 - g) / 4.0
    else:
        gpp = (1.0 - g) - (1.0/g) * gp**2 - (2.0/r) * gp
    return [gp, gpp]

def solve_ode(g0, r_max=400.0, n_pts=160000):
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=1e-13, atol=1e-15, max_step=0.015)
    return sol

def extract_atail(sol, r_min=100.0, r_max=350.0):
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    return math.sqrt(coef[0]**2 + coef[1]**2)

def measure_eta(delta):
    sol = solve_ode(1.0 + delta)
    if not sol.success: return None
    A = extract_atail(sol)
    return A / abs(delta)

# ==================== HIGH-PRECISION MEASUREMENT ====================
print("=" * 78)
print("  IDENTYFIKACJA CLOSED-FORM dla alpha_3 (O(delta^2) w eta)")
print("=" * 78)

C1_THEORY = 1.0 - math.log(3.0)/4.0
print(f"\n  c1 = 1 - ln(3)/4 = {C1_THEORY:.14f}  (z breakthrough)")
print(f"  alpha_2 = c1/2   = {C1_THEORY/2:.14f}")

# Gesty grid malych delta dla precyzyjnej ekstrapolacji
deltas = np.array([0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.012, 0.02, 0.03, 0.05])

print(f"\n  Pomiar eta(+delta) i eta(-delta):")
print(f"  {'delta':>8s}  {'eta_exc':>16s}  {'eta_def':>16s}  {'sym-1':>14s}")

meas = []
for d in deltas:
    e_p = measure_eta(+d)
    e_m = measure_eta(-d)
    sym = (e_p + e_m) / 2.0
    meas.append((d, e_p, e_m, sym))
    print(f"  {d:8.5f}  {e_p:16.12f}  {e_m:16.12f}  {sym-1:14.6e}")

deltas_arr = np.array([m[0] for m in meas])
sym_arr = np.array([m[3] for m in meas])

# ==================== RICHARDSON EXTRAPOLATION ====================
# sym(delta) = 1 + alpha_3*delta^2 + alpha_5*delta^4 + alpha_7*delta^6 + ...
# Uzyj najmniejszych delta, fit polynomial w delta^2

print(f"\n{'-'*78}")
print(f"  Polynomial fit: sym(delta) - 1 = a3*d^2 + a5*d^4 + a7*d^6")
print(f"{'-'*78}")

y = sym_arr - 1.0

# Multiple fit strategies
# (1) All points, cubic in d^2
X1 = np.column_stack([deltas_arr**2, deltas_arr**4, deltas_arr**6])
c1_all = np.linalg.lstsq(X1, y, rcond=None)[0]
print(f"\n  Fit all {len(deltas)} punktow, poly w d^2 stopnia 3:")
print(f"    alpha_3 = {c1_all[0]:.12f}")
print(f"    alpha_5 = {c1_all[1]:+.12f}")
print(f"    alpha_7 = {c1_all[2]:+.12f}")

# (2) Smallest 6 points, quadratic in d^2
idx = np.argsort(deltas_arr)[:6]
Xs = np.column_stack([deltas_arr[idx]**2, deltas_arr[idx]**4])
c1_s = np.linalg.lstsq(Xs, y[idx], rcond=None)[0]
print(f"\n  Fit 6 smallest, poly w d^2 stopnia 2:")
print(f"    alpha_3 = {c1_s[0]:.12f}")
print(f"    alpha_5 = {c1_s[1]:+.12f}")

# (3) Smallest 4, linear in d^2
idx2 = np.argsort(deltas_arr)[:4]
Xl = np.column_stack([deltas_arr[idx2]**2])
c1_l = np.linalg.lstsq(Xl, y[idx2], rcond=None)[0]
print(f"\n  Fit 4 smallest, linear w d^2:")
print(f"    alpha_3 = {c1_l[0]:.12f}")

# (4) Richardson: c1(d1)*d2^2 - c1(d2)*d1^2 / (d2^2 - d1^2)
# 2-point Richardson dla sym/d^2
r_est = y / deltas_arr**2
print(f"\n  Richardson: sym/d^2 dla roznych d:")
for i, d in enumerate(deltas_arr):
    print(f"    d={d:7.5f}: sym/d^2 = {r_est[i]:.12f}")

# Fit r_est = alpha_3 + alpha_5 * d^2 + ...
# Use small d only
idx_r = np.argsort(deltas_arr)[:5]
X_r = np.column_stack([np.ones(5), deltas_arr[idx_r]**2])
c_r = np.linalg.lstsq(X_r, r_est[idx_r], rcond=None)[0]
alpha3_rich = c_r[0]
print(f"\n  Richardson extrapolation alpha_3 (d->0) = {alpha3_rich:.12f}")

# ==================== COMBINATORIAL SCAN ====================
# Skan: alpha_3 = a*pi^2 + b*ln(3) + c*ln(2) + d*(1-c1) + e  dla rationals a,b,c,d,e
print(f"\n{'-'*78}")
print(f"  COMBINATORIAL SCAN: alpha_3 jako kombinacja ladnych stalych")
print(f"{'-'*78}")

alpha3_target = alpha3_rich  # Richardson-extrapolated value

# Kandydaci: stale podstawowe
BASES = {
    '1': 1.0,
    'pi': math.pi,
    'pi^2': math.pi**2,
    'pi^3': math.pi**3,
    'ln(2)': math.log(2),
    'ln(3)': math.log(3),
    'ln(5)': math.log(5),
    '(ln3)^2': math.log(3)**2,
    'c1=1-ln3/4': C1_THEORY,
    '1-c1=ln3/4': 1-C1_THEORY,
    'c1^2': C1_THEORY**2,
    '(1-c1)^2': (1-C1_THEORY)**2,
    'c1*(1-c1)': C1_THEORY*(1-C1_THEORY),
    'zeta(3)': 1.2020569031595942,
    'Catalan': 0.9159655941772190,
    'gamma': 0.5772156649015329,  # Euler-Mascheroni
}

# 1-term match: alpha_3 = const/N for small N
print(f"\n  1-term: alpha_3 = const/N (N in 1..200)")
best_1 = []
for name, val in BASES.items():
    for N in range(1, 201):
        est = val / N
        rel = abs(est - alpha3_target) / abs(alpha3_target)
        if rel < 1e-3:
            best_1.append((name, N, est, rel))
for name, N, est, rel in sorted(best_1, key=lambda x: x[3])[:20]:
    print(f"    {name:15s} / {N:4d} = {est:.12f}  rel={rel:.3e}")

# 2-term match: alpha_3 = a*const1/N1 + b*const2/N2 (a,b in {-1,1,2,-2,1/2,-1/2})
print(f"\n  2-term: alpha_3 = a*B1/N1 + b*B2/N2 (a,b in simple rationals)")
print(f"  ({'pair':45s}  value              rel_err)")

best_2 = []
simple_rats = [1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 3.0, 1/3, -1/3]
base_list = list(BASES.items())
for (n1, v1), (n2, v2) in [(a, b) for i, a in enumerate(base_list) for j, b in enumerate(base_list) if i < j]:
    for N1 in [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48]:
        for N2 in [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48]:
            for a, b in product(simple_rats, simple_rats):
                est = a * v1 / N1 + b * v2 / N2
                if abs(est) > 10: continue  # avoid huge cancellations
                rel = abs(est - alpha3_target) / abs(alpha3_target)
                if rel < 1e-4:
                    pair_str = f"{a:+.2g}*{n1}/{N1} + {b:+.2g}*{n2}/{N2}"
                    best_2.append((pair_str, est, rel))

best_2 = sorted(set(best_2), key=lambda x: x[2])[:15]
for pair_str, est, rel in best_2:
    print(f"  {pair_str:45s}  {est:.12f}  {rel:.3e}")

# Special candidates exploiting c1 structure
print(f"\n  Hipotezy oparte na c1 i ln(3):")
special = {
    'c1^2 / 8':                C1_THEORY**2 / 8,
    'c1^2 / 6':                C1_THEORY**2 / 6,
    'c1/8':                    C1_THEORY/8,
    '(1-c1)^2':                (1-C1_THEORY)**2,
    'c1*(1-c1)':               C1_THEORY*(1-C1_THEORY),
    'c1*(1-c1)/2':             C1_THEORY*(1-C1_THEORY)/2,
    'c1^3':                    C1_THEORY**3,
    'ln(3)/12':                math.log(3)/12,
    'ln(3)^2/16':              math.log(3)**2/16,
    'ln(3)^2/12':              math.log(3)**2/12,
    '(1-c1)*ln(3)':            (1-C1_THEORY)*math.log(3),
    'c1*ln(3)/9':              C1_THEORY*math.log(3)/9,
    'c1*ln(3)/8':              C1_THEORY*math.log(3)/8,
    'ln(3)*c1^2/3':            math.log(3)*C1_THEORY**2/3,
    '(ln 3 * pi^2) / 120':     math.log(3)*math.pi**2/120,
    '(ln 3)^3 / 16':           math.log(3)**3/16,
    '1/12 + 1/120':            1/12 + 1/120,
    '1/12 + c1^2/100':         1/12 + C1_THEORY**2/100,
    '2/9 - c1/6':              2/9 - C1_THEORY/6,
    'pi^2/110':                math.pi**2/110,
    'pi^2/11/10':              math.pi**2/(11*10),
    '(1 - 2*ln(3)/3)/8':       (1 - 2*math.log(3)/3)/8,
    '(2 - ln 3)/10':           (2 - math.log(3))/10,
    '(3/8)*(1 - 2*ln(3)/3)':   (3.0/8.0)*(1 - 2*math.log(3)/3),
    '1/12 + (ln 3)^2 / 60':    1/12 + math.log(3)**2/60,
    'pi^2/64 - ln(3)/16':      math.pi**2/64 - math.log(3)/16,
    'c1*(1-c1)/8':             C1_THEORY*(1-C1_THEORY)/8,
    '(2-c1)/16':               (2-C1_THEORY)/16,
    '3*c1/16 - 1/32':          3*C1_THEORY/16 - 1/32,
    '(1 - c1/2)/16':           (1 - C1_THEORY/2)/16,
    '(3/2 - c1)^2 / 8':        (1.5 - C1_THEORY)**2 / 8,
}

print(f"\n  {'Kandydat':>30s}  {'Wartosc':>16s}  {'Diff':>14s}  {'Rel':>10s}")
print(f"  {'-'*30}  {'-'*16}  {'-'*14}  {'-'*10}")
sorted_s = sorted(special.items(), key=lambda kv: abs(kv[1]-alpha3_target))
for name, val in sorted_s[:20]:
    d = val - alpha3_target
    rel = abs(d) / abs(alpha3_target)
    print(f"  {name:>30s}  {val:16.12f}  {d:+14.4e}  {rel:10.3e}")

print(f"\n{'-'*78}")
print(f"  PODSUMOWANIE")
print(f"{'-'*78}")
print(f"  alpha_3 (Richardson, d->0)     = {alpha3_rich:.12f}")
print(f"  alpha_3 (wszystkie punkty)     = {c1_all[0]:.12f}")
print(f"  alpha_3 (6 najmniejszych)      = {c1_s[0]:.12f}")
print(f"  alpha_3 (4 najmniejszych, lin) = {c1_l[0]:.12f}")

spread = max(c1_all[0], c1_s[0], c1_l[0], alpha3_rich) - min(c1_all[0], c1_s[0], c1_l[0], alpha3_rich)
print(f"\n  Spread miedzy metodami: {spread:.4e}")
print(f"  (ten poziom = systematyczny bias pomiaru ODE przy r_max=400)")
