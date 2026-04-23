#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c2_alpha5_identify.py
=========================

Identyfikacja alpha_5 (coefficient przy delta^4 w eta_sym) po sukcesie z alpha_3.

Struktura znana:
  eta_sym(delta) = (eta_exc + eta_def)/2 = 1 + alpha_3*delta^2 + alpha_5*delta^4 + ...
  alpha_2 = c1/2 = (1 - ln(3)/4)/2          [DOKLADNIE, analitycznie]
  alpha_3 = pi^2/110 = pi^2/128 + 9*pi^2/7040 [IDENTIFIED + NUMERIC VERIFIED]
  alpha_5 = ???                                [CEL]

Strategia:
  - Wysoka precyzja (r_max=1600) dla wielu delta
  - Fit do delta^6
  - Pokusic sie o alpha_5 (a moze tez alpha_7)
  - Test kandydatow: pi^2, pi^4, ln(3), ...

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

def solve_ode(g0, r_max=1600.0, n_pts=640000):
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=1e-13, atol=1e-15, max_step=0.015)
    return sol

def extract_atail(sol, r_min=200.0, r_max=1550.0):
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

# Mnie deltas dla szerokiego fitu
deltas = np.array([0.003, 0.005, 0.008, 0.012, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1])

alpha3_target = math.pi**2 / 110.0
c1_th = 1.0 - math.log(3.0)/4.0

print("=" * 78)
print("  alpha_5 IDENTIFIKACJA (po alpha_3 = pi^2/110)")
print("=" * 78)

sym_arr = []
print(f"\n  Pomiar sym(delta) - 1 (r_max=1600, n_pts=640k):")
print(f"  {'delta':>8s}  {'sym-1':>14s}  {'/delta^2':>14s}")
for d in deltas:
    e_p = measure_eta(+d)
    e_m = measure_eta(-d)
    sym = (e_p + e_m) / 2.0
    sym_arr.append(sym - 1.0)
    print(f"  {d:8.5f}  {sym-1:14.10e}  {(sym-1)/d**2:14.10f}")

y = np.array(sym_arr)

# Fit y = a3*d^2 + a5*d^4 + a7*d^6
print(f"\n  FIT POLYNOMIAL (d^2, d^4, d^6):")
X = np.column_stack([deltas**2, deltas**4, deltas**6])
coef, *_ = np.linalg.lstsq(X, y, rcond=None)
alpha_3 = coef[0]
alpha_5 = coef[1]
alpha_7 = coef[2]
print(f"    alpha_3 = {alpha_3:.12f}")
print(f"    alpha_5 = {alpha_5:.12f}")
print(f"    alpha_7 = {alpha_7:.12f}")
print(f"    (target alpha_3 = pi^2/110 = {alpha3_target:.12f}, diff = {alpha_3-alpha3_target:+.3e})")

# Alternative: fit with CONSTRAINED alpha_3 = pi^2/110
# y - alpha3_target*d^2 = alpha_5*d^4 + alpha_7*d^6
print(f"\n  FIT Z alpha_3 = pi^2/110 ZADANE:")
y_corr = y - alpha3_target * deltas**2
X2 = np.column_stack([deltas**4, deltas**6])
coef2, *_ = np.linalg.lstsq(X2, y_corr, rcond=None)
alpha_5_c = coef2[0]
alpha_7_c = coef2[1]
print(f"    alpha_5 = {alpha_5_c:.12f}")
print(f"    alpha_7 = {alpha_7_c:.12f}")

# Residuals dla bootstrap uncertainty
pred = X2 @ coef2
resid = y_corr - pred
rms = np.sqrt(np.mean(resid**2))
print(f"    RMS residual = {rms:.3e}")

# Kandydaci closed form dla alpha_5
print(f"\n{'='*78}")
print(f"  KANDYDACI CLOSED FORM dla alpha_5")
print(f"{'='*78}")
print(f"\n  alpha_5 (constrained fit) = {alpha_5_c:.10f}")
print(f"  alpha_5 (unconstrained)   = {alpha_5:.10f}")
print(f"  (roznica = {alpha_5 - alpha_5_c:+.3e})")

target_a5 = alpha_5_c

pi = math.pi
ln3 = math.log(3.0)
ln2 = math.log(2.0)

cands = {}

# Naturalne wzorce: rational * pi^k lub wyrazone przez I_sin, I_cos
# From the structure of pertubation, expect:
# alpha_5 = some combination of I_sin^2, I_cos^2, sort of products
# I_cos = 1/2 - ln(3)/8, I_sin = -pi/8

# Try rationals of pi, pi^2, pi^3, pi^4
for n in [pi, pi**2, pi**3, pi**4, ln3, ln3**2, ln3**3]:
    for den in range(10, 5000, 1):
        v = n/den
        if abs(v - target_a5) < 1e-4:
            cands[f'{n:.6g}/{den}'] = v

# Explicit structural candidates based on expansion
# Leading O(delta^4) from eta expansion:
# (1+x)^0.5 expansion: if A^2/delta^2 = 1 + 2·I_cos·delta + (I_cos^2+I_sin^2)delta^2 + ...
# Then (A/delta)^2 = 1 + 2·a2·d + 2·a3·d^2 + 2·a5·d^4 + (a2^2 + 2·a3) d^2 + ... hmm

# Let's use cleaner structure:
# eta_sym(d) = 1 + a3*d^2 + a5*d^4 + ...
# Physical prediction from pertubation: a5 = f(I_cos, I_sin, P_cos, P_sin, Q_cos, ...)

# Pure candidate sweep
explicit_cands = {
    'pi^4/1024':        pi**4/1024,
    'pi^4/2048':        pi**4/2048,
    'pi^4/3430':        pi**4/3430,
    'pi^2/352':         pi**2/352,
    'pi^4/3600':        pi**4/3600,
    'pi^2/500':         pi**2/500,
    'pi^2/400':         pi**2/400,
    'pi^2/512':         pi**2/512,
    '3*pi^2/512':       3*pi**2/512,
    'pi/128':           pi/128,
    'I_sin^2 * I_cos':  (-pi/8)**2 * (0.5 - ln3/8),
    'I_cos^2 * pi':     (0.5 - ln3/8)**2 * pi,
    'I_cos^4 / 2':      (0.5 - ln3/8)**4 / 2,
    '(1 - ln3/4)^2 /4': (1-ln3/4)**2 / 4,
    'c1^2 / 30':        c1_th**2 / 30,
    'c1^2 / 40':        c1_th**2 / 40,
    'c1^4 / 8':         c1_th**4 / 8,
    'c1 * pi / 50':     c1_th * pi / 50,
    'pi^2 / (110 * 3)': pi**2 / 330,
    'pi^2 / (110 * 4)': pi**2 / 440,
    'pi^4 / (110^2)':   pi**4 / 12100,
    '(pi^2/110)^2 /2':  (pi**2/110)**2 / 2,
}

print(f"\n  Strukturalne kandydaty:")
print(f"  {'Kandydat':30s}  {'Wartosc':>14s}  {'Diff':>12s}  {'Rel':>10s}")
print(f"  {'-'*30}  {'-'*14}  {'-'*12}  {'-'*10}")
for name, v in sorted(explicit_cands.items(), key=lambda kv: abs(kv[1]-target_a5)):
    d = v - target_a5
    rel = abs(d) / abs(target_a5) if target_a5 != 0 else float('inf')
    marker = " <--" if abs(d) < 1e-3 else ""
    print(f"  {name:30s}  {v:14.10f}  {d:+12.3e}  {rel:10.2e}{marker}")
