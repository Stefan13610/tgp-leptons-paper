#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c2_perturbative_coeffs.py
=============================

Kontynuacja breakthrough c1 = 1 - ln(3)/4: obliczamy WYZSZE wspolczynniki
rozwiniecia perturbacyjnego A_tail(1+delta) dla malych delta.

TEORIA:
  g(r; g0) = 1 + delta*f(r) + delta^2*h(r) + delta^3*p(r) + ...
  gdzie delta = g0 - 1, a f, h, p sa rzadami rozwiniecia.

  A_tail(1+delta) = delta * [1 + alpha_2*delta + alpha_3*delta^2 + alpha_4*delta^3 + ...]

  eta(delta) = A_tail(1+delta) / delta  (signed, for excess)
           = 1 + alpha_2*delta + alpha_3*delta^2 + alpha_4*delta^3 + ...

  Dla deficit (delta -> -delta, ale eta liczone jako A_tail/|delta|):
  eta_def(delta) = A_tail(1-delta) / delta  (delta > 0)
                = 1 - alpha_2*delta + alpha_3*delta^2 - alpha_4*delta^3 + ...

  ASYM: (eta_exc - eta_def)/delta = 2*alpha_2 + 2*alpha_4*delta^2 + O(delta^4)
  SYM:  (eta_exc + eta_def)/2    = 1 + alpha_3*delta^2 + alpha_5*delta^4 + O(delta^6)

  Z breakthrough: c1 == 2*alpha_2 = 1 - ln(3)/4  (UDOWODNIONE)

CEL:
  Zmierzyc numerycznie:
    - alpha_2 (weryfikacja c1 = 1 - ln(3)/4)
    - alpha_3 (NOWA stala, do identyfikacji closed form)
    - alpha_4 (korekcja do c1)

  Sprobowac zidentyfikowac alpha_3 w formie zamknietej (pi^2, ln(2), 1/12, ...).

HIPOTEZY KANDYDATSKIE dla alpha_3:
  - pi^2/48  ≈ 0.2056
  - 1/6      ≈ 0.1667
  - ln(2)/4  ≈ 0.1733
  - 1/(4π)   ≈ 0.0796
  - (ln(3))^2 / 16 ≈ 0.0754
  - 7/60     ≈ 0.1167

Author: Claudian
Date: 2026-04-16
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp

# ==================== ODE solver ====================
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

def extract_atail_signed(sol, r_min=100.0, r_max=350.0, sign_delta=+1):
    """
    Wyciag A_tail JAKO LICZBE ZE ZNAKIEM.

    Obserwacja: dla deficit (delta<0) i excess (delta>0), tail u(r) = (g-1)*r
    zaczyna sie z przeciwnymi znakami (u blisko r_match).
    Przy dopasowaniu A*cos(r-phi) normalizujemy tak, ze A_tail = sign(delta) * sqrt(B^2+C^2).

    Wrodomosc: A_tail(1+delta) jest nieparzysta w leading order (delta wiodacy).
    """
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    B, C = coef
    A_mag = math.sqrt(B*B + C*C)
    # Sign: dominant sign of u_f (mean over window)
    sign = 1.0 if np.mean(u_f) * sign_delta > 0 else (-1.0 if np.mean(u_f) * sign_delta < 0 else 1.0)
    # Lepsze: sign z pierwszej oscylacji (r blisko r_min)
    # Ale dla leading-order pewniejsze: zalozmy A_tail \propto sign_delta
    # aby eta(delta) = A_tail / delta byla CIAGLA w delta -> 0
    # Zwrocmy MAGNITUDE; znak zarzadzany przez callera
    y_hat = B * np.cos(r_f) + C * np.sin(r_f)
    rmse = float(np.sqrt(np.mean((u_f - y_hat)**2)))
    return A_mag, rmse / max(A_mag, 1e-15)

def get_eta(delta):
    """eta(delta) = A_tail(1+delta) / |delta| (dodatnie; dla eta_def i eta_exc)"""
    g0 = 1.0 + delta
    sol = solve_ode(g0)
    if not sol.success:
        return None, None
    A, q = extract_atail_signed(sol, sign_delta=np.sign(delta))
    if q > 1e-3:
        return None, q
    return A / abs(delta), q

# ==================== MEASUREMENT ====================
print("=" * 76)
print("  R6 Sciezka D: c2 i wyzsze - wspolczynniki rozwiniecia eta(delta)")
print("=" * 76)

C1_THEORY = 1.0 - math.log(3.0)/4.0
ALPHA2_THEORY = C1_THEORY / 2.0

print(f"\n  Z breakthrough 2026-04-16:")
print(f"    c1 = 2*alpha_2 = 1 - ln(3)/4 = {C1_THEORY:.12f}")
print(f"    alpha_2        = {ALPHA2_THEORY:.12f}")

# Measurement grid: small delta, both excess and deficit
deltas = np.array([0.002, 0.003, 0.005, 0.008, 0.012, 0.018, 0.025, 0.035, 0.050])

print(f"\n{'-'*76}")
print(f"  Pomiar A_tail(1+/-delta) z rtol=1e-13, n_pts=160k:")
print(f"{'-'*76}")
print(f"\n  {'delta':>8s}  {'eta_exc':>16s}  {'eta_def':>16s}  {'sym':>16s}  {'asym/d':>16s}  {'q':>8s}")
print(f"  {'-'*8}  {'-'*16}  {'-'*16}  {'-'*16}  {'-'*16}  {'-'*8}")

measurements = []
for d in deltas:
    eta_e, q_e = get_eta(+d)
    eta_d, q_d = get_eta(-d)
    if eta_e is None or eta_d is None:
        print(f"  {d:8.5f}  FAILED")
        continue
    sym = (eta_e + eta_d) / 2.0
    asym_over_d = (eta_e - eta_d) / d
    q = max(q_e, q_d)
    measurements.append((d, eta_e, eta_d, sym, asym_over_d))
    print(f"  {d:8.5f}  {eta_e:16.12f}  {eta_d:16.12f}  {sym:16.12f}  {asym_over_d:16.12f}  {q:8.2e}")

# ==================== FIT POLYNOMIALS ====================
print(f"\n{'-'*76}")
print(f"  Fit: eta_exc(delta) = 1 + alpha_2*delta + alpha_3*delta^2 + alpha_4*delta^3 + alpha_5*delta^4")
print(f"{'-'*76}")

deltas_arr = np.array([m[0] for m in measurements])
eta_exc_arr = np.array([m[1] for m in measurements])
eta_def_arr = np.array([m[2] for m in measurements])
sym_arr = np.array([m[3] for m in measurements])
asym_arr = np.array([m[4] for m in measurements])

# Fit eta_exc - 1 = a2*d + a3*d^2 + a4*d^3 + a5*d^4 (order 4)
# Uzyc najmniejszych delta do minimalizacji wyzszych rzedow
y_exc = eta_exc_arr - 1.0
y_def = eta_def_arr - 1.0

# Build Vandermonde
def poly_fit_coeffs(deltas, y_values, max_order=4):
    """Fit y = c1*d + c2*d^2 + ... + c_max*d^max, return coeffs in order [c1, c2, ...]"""
    X = np.column_stack([deltas**k for k in range(1, max_order+1)])
    coef, *_ = np.linalg.lstsq(X, y_values, rcond=None)
    return coef

# Fit ze wszystkich deltas
coef_exc = poly_fit_coeffs(deltas_arr, y_exc, max_order=4)
coef_def = poly_fit_coeffs(deltas_arr, y_def, max_order=4)

print(f"\n  eta_exc(delta) - 1 coefficients:")
for i, c in enumerate(coef_exc, start=1):
    print(f"    alpha_{i+1} = {c:+.10f}")

print(f"\n  eta_def(delta) - 1 coefficients (for delta > 0):")
for i, c in enumerate(coef_def, start=1):
    print(f"    c_def_{i} = {c:+.10f}")

# Extract alpha_2 (z excess: wspolczynnik przy delta^1)
alpha2_exc = coef_exc[0]
alpha3_exc = coef_exc[1]
alpha4_exc = coef_exc[2]
alpha5_exc = coef_exc[3]

# alpha_2 z deficit powinno byc -alpha_2_exc (odwrotny znak)
alpha2_def = coef_def[0]

print(f"\n{'-'*76}")
print(f"  alpha_2 weryfikacja: eta_exc daje alpha_2 = {alpha2_exc:.10f}")
print(f"                       eta_def daje -alpha_2 = {alpha2_def:.10f}")
print(f"                       teoria:     alpha_2  = {ALPHA2_THEORY:.10f}")
print(f"  diff (exc)         = {alpha2_exc - ALPHA2_THEORY:+.4e}")
print(f"  diff (def, neg)    = {-alpha2_def - ALPHA2_THEORY:+.4e}")

# ==================== ALPHA_3 IDENTIFICATION ====================
print(f"\n{'-'*76}")
print(f"  alpha_3 (z excess i z sym): identyfikacja closed form")
print(f"{'-'*76}")

# alpha_3 z sym(delta) = 1 + alpha_3*delta^2 + alpha_5*delta^4 + ...
sym_y = sym_arr - 1.0
X_sym = np.column_stack([deltas_arr**2, deltas_arr**4])
coef_sym, *_ = np.linalg.lstsq(X_sym, sym_y, rcond=None)
alpha3_sym = coef_sym[0]
alpha5_sym = coef_sym[1]

# alpha_3 z eta_exc
print(f"\n  alpha_3 (z eta_exc poly fit)  = {alpha3_exc:.10f}")
print(f"  alpha_3 (z sym = (exc+def)/2) = {alpha3_sym:.10f}")
print(f"  alpha_5 (z sym)               = {alpha5_sym:.10f}")

alpha3_measured = alpha3_sym  # use symmetric - more robust

# Scan closed-form candidates
print(f"\n  Skan kandydatow closed-form dla alpha_3:")
print(f"  {'Kandydat':>30s}  {'Wartosc':>14s}  {'Diff od alpha_3':>16s}  {'Rel':>10s}")
print(f"  {'-'*30}  {'-'*14}  {'-'*16}  {'-'*10}")

candidates_3 = {
    "pi^2 / 48":                math.pi**2 / 48,
    "pi^2 / 64":                math.pi**2 / 64,
    "pi^2 / 96":                math.pi**2 / 96,
    "pi^2 / 120":               math.pi**2 / 120,
    "1/6":                      1/6,
    "1/12":                     1/12,
    "1/8":                      1/8,
    "1/16":                     1/16,
    "1/32":                     1/32,
    "ln(2)/4":                  math.log(2)/4,
    "ln(2)/8":                  math.log(2)/8,
    "ln(2)/12":                 math.log(2)/12,
    "ln(3)^2 / 16":             math.log(3)**2 / 16,
    "ln(3)^2 / 32":             math.log(3)**2 / 32,
    "ln(3)/16":                 math.log(3)/16,
    "ln(3)/24":                 math.log(3)/24,
    "1/(4*pi)":                 1/(4*math.pi),
    "1/(8*pi)":                 1/(8*math.pi),
    "c1^2 / 2":                 C1_THEORY**2 / 2,
    "c1^2 / 4":                 C1_THEORY**2 / 4,
    "c1^2 / 8":                 C1_THEORY**2 / 8,
    "(1-c1)^2":                 (1-C1_THEORY)**2,
    "(1-c1)":                   (1-C1_THEORY),
    "(1-c1)/2":                 (1-C1_THEORY)/2,
    "(1-c1)/4":                 (1-C1_THEORY)/4,
    "1 - c1":                   1 - C1_THEORY,  # = ln(3)/4 = 0.2747
    "ln(5)/16":                 math.log(5)/16,
    "pi/32":                    math.pi/32,
    "pi/48":                    math.pi/48,
    "7/60":                     7/60,
    "7/72":                     7/72,
    "3/32":                     3/32,
    "11/96":                    11/96,
    "(3-2*ln(3))/8":            (3 - 2*math.log(3))/8,
    "(1 - ln(3)/3)/4":          (1 - math.log(3)/3)/4,
    "ln(3)^2 / 48":             math.log(3)**2 / 48,
}

sorted_cands = sorted(candidates_3.items(), key=lambda kv: abs(kv[1] - alpha3_measured))
for name, val in sorted_cands[:15]:
    d = val - alpha3_measured
    rel = abs(d) / max(abs(alpha3_measured), 1e-12)
    print(f"  {name:>30s}  {val:14.10f}  {d:+16.4e}  {rel:10.2e}")

# ==================== REFINED FIT (smallest deltas only) ====================
print(f"\n{'-'*76}")
print(f"  Refined fit: uzywamy tylko 5 najmniejszych delta dla wiekszej precyzji")
print(f"{'-'*76}")

# Top 5 smallest deltas
N_small = 5
idx_sorted = np.argsort(deltas_arr)[:N_small]
deltas_small = deltas_arr[idx_sorted]
sym_small = sym_arr[idx_sorted] - 1.0

X_ref = np.column_stack([deltas_small**2, deltas_small**4])
coef_ref, *_ = np.linalg.lstsq(X_ref, sym_small, rcond=None)
alpha3_refined = coef_ref[0]
alpha5_refined = coef_ref[1]
print(f"\n  alpha_3 (refined)  = {alpha3_refined:.12f}")
print(f"  alpha_5 (refined)  = {alpha5_refined:+.12f}")

# Same dla ASYM (-> c1 i korekcja)
asym_small = asym_arr[idx_sorted]
coef_asym, *_ = np.linalg.lstsq(X_ref, asym_small, rcond=None)
# asym/d = 2*alpha_2 + 2*alpha_4*delta^2 + ...
c1_refined = coef_asym[0]         # to jest 2*alpha_2
correction1 = coef_asym[1]        # to jest 2*alpha_4
# Ale wait: we fit X_ref = [d^2, d^4] but asym/d = 2a2 + 2a4*d^2 + 2a6*d^4 - poly w d^2
# Zreinteresowac:
X_asym = np.column_stack([np.ones_like(deltas_small), deltas_small**2])
coef_asym2, *_ = np.linalg.lstsq(X_asym, asym_small, rcond=None)
c1_refined_v2 = coef_asym2[0]
print(f"\n  c1 (refined from asym/d -> delta=0) = {c1_refined_v2:.12f}")
print(f"  Target: 1 - ln(3)/4                = {C1_THEORY:.12f}")
print(f"  Diff                               = {c1_refined_v2 - C1_THEORY:+.4e}")

# alpha_3 refined match
print(f"\n  Skan dla alpha_3 (refined) = {alpha3_refined:.10f}:")
print(f"  {'Kandydat':>30s}  {'Wartosc':>14s}  {'Rel err':>10s}")
for name, val in sorted_cands[:6]:
    rel = abs(val - alpha3_refined) / max(abs(alpha3_refined), 1e-12)
    print(f"  {name:>30s}  {val:14.10f}  {rel:10.2e}")

print(f"\n{'-'*76}")
print(f"  PODSUMOWANIE:")
print(f"{'-'*76}")
print(f"  alpha_2 (c1/2) zmierzone  = {alpha2_exc:.10f}")
print(f"  alpha_2 teoria            = {ALPHA2_THEORY:.10f}")
print(f"  Agreement                 = {abs(alpha2_exc-ALPHA2_THEORY):.1e}")
print()
print(f"  alpha_3 zmierzone (sym)   = {alpha3_refined:.10f}")
print(f"  Najblizszy kandydat       = {sorted_cands[0][0]} = {sorted_cands[0][1]:.10f}")
print(f"  Diff                      = {sorted_cands[0][1] - alpha3_refined:+.4e}")
