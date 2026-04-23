#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c1_high_precision.py
=======================

Test wysoko-precyzyjny hipotezy:

    c1 = 1 - ln(3)/4   <=>   exp(4*(1-c1)) = 3

gdzie c1 to STALA asymetrii deficit/excess w ogonie ODE solitonowego TGP:

    eta(delta) = A_tail(1 + delta) / |delta|
    (eta_exc(delta) - eta_def(delta)) / delta -> c1   dla  delta -> 0

Z dotychczasowych pomiarow (r6_eta_koide_attack.py, precyzja ~10^-5):
    c1 ≈ 0.72538

Kandydat zamkniety:
    1 - ln(3)/4 = 0.72534693  (diff = +3.3e-5 od pomiaru)

Czy to przypadek, czy relacja scisla?

Metoda:
  1. Rozwiaz ODE z bardzo wysoka precyzja (rtol=1e-13)
  2. Uzyj DOLPHIN-styl extrapolacji Richardson: probki c1(delta) dla
     delta = 0.01, 0.005, 0.002, 0.001, i extrapoluj do delta=0
  3. Porownaj z 1 - ln(3)/4 na poziomie 10^-7

Physical interpretation (jesli PASS):
  (1-c1) = ln(3)/4 wiaze asymetrie deficit/excess ODE z N=3 generacjami
  poprzez ln(3) - charakterystyka 3-poziomowej entropii (Shannon H = ln(3)
  dla uniform distribution na 3 stanach).

Author: Claudian
Date: 2026-04-16
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math
from scipy.integrate import solve_ivp

# Target closed-form candidate
C1_HYPOTHESIS = 1.0 - math.log(3.0) / 4.0   # = 0.7253469278...

print("=" * 72)
print("  R6.8: Wysoko-precyzyjny test hipotezy c1 = 1 - ln(3)/4")
print("=" * 72)
print(f"\n  Hipoteza:  c1 = 1 - ln(3)/4 = {C1_HYPOTHESIS:.12f}")
print(f"             exp(4(1-c1)) = 3 exactly")
print(f"             (interpretacja: log-3 = Shannon entropy dla N=3 generacji)")

# ------------------------------------------------------------------ #
# High-precision ODE solver                                           #
# ------------------------------------------------------------------ #

def rhs(r, y):
    g, gp = y
    if g < 1e-12:
        g = 1e-12
    if r < 1e-13:
        gpp = (1.0 - g) / 4.0
    else:
        gpp = (1.0 - g) - (1.0 / g) * gp**2 - (2.0 / r) * gp
    return [gp, gpp]


def solve_ode(g0, r_max=400.0, n_pts=120000):
    r_eval = np.linspace(1e-10, r_max, n_pts)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853',           # 8th-order Runge-Kutta
                    t_eval=r_eval,
                    rtol=1e-13, atol=1e-15,
                    max_step=0.02)
    return sol


def extract_atail(sol, r_min=100.0, r_max=350.0):
    """Fit u(r) = (g-1)*r = B*cos(r) + C*sin(r); A_tail = sqrt(B^2+C^2)"""
    r = sol.t
    g = sol.y[0]
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, *_ = np.linalg.lstsq(X, u_f, rcond=None)
    A = math.sqrt(coef[0]**2 + coef[1]**2)
    # Residual quality check
    y_hat = coef[0] * np.cos(r_f) + coef[1] * np.sin(r_f)
    rmse = float(np.sqrt(np.mean((u_f - y_hat)**2)))
    return A, rmse / max(A, 1e-12)


def get_atail(g0):
    sol = solve_ode(g0)
    if not sol.success:
        return None
    A, q = extract_atail(sol)
    if q > 1e-3:
        # degraded fit; expand range
        A, q = extract_atail(sol, r_min=80.0, r_max=380.0)
    return A


# ------------------------------------------------------------------ #
# Compute c1(delta) for small delta                                   #
# ------------------------------------------------------------------ #

print(f"\n{'-'*72}")
print("  c1(delta) = (eta_exc - eta_def) / delta   (-> c1 gdy delta -> 0)")
print(f"{'-'*72}")

# Uzywamy kilku malych delta aby Richardson-extrapolowac do delta=0
deltas = [0.05, 0.03, 0.02, 0.015, 0.01, 0.008, 0.006, 0.004, 0.002]

print(f"\n  {'delta':>8s}  {'eta_def':>14s}  {'eta_exc':>14s}  {'c1(delta)':>14s}  {'q_def':>8s}  {'q_exc':>8s}")
print(f"  {'-'*8}  {'-'*14}  {'-'*14}  {'-'*14}  {'-'*8}  {'-'*8}")

results = []
for d in deltas:
    A_def = get_atail(1.0 - d)
    A_exc = get_atail(1.0 + d)
    if A_def is None or A_exc is None:
        print(f"  {d:8.5f}   FAILED")
        continue
    eta_def = A_def / d
    eta_exc = A_exc / d
    c1_local = (eta_exc - eta_def) / d
    results.append((d, c1_local, eta_def, eta_exc))
    # Szybki quality
    sol_def = solve_ode(1.0 - d)
    _, q_def = extract_atail(sol_def)
    sol_exc = solve_ode(1.0 + d)
    _, q_exc = extract_atail(sol_exc)
    print(f"  {d:8.5f}  {eta_def:14.10f}  {eta_exc:14.10f}  "
          f"{c1_local:14.10f}  {q_def:8.2e}  {q_exc:8.2e}")

# ------------------------------------------------------------------ #
# Richardson extrapolation: c1(delta) = c1_0 + c2 * delta^2 + ...      #
# ------------------------------------------------------------------ #

print(f"\n{'-'*72}")
print("  Richardson extrapolacja do delta = 0")
print(f"{'-'*72}")

results.sort(key=lambda x: x[0])  # Sort by delta ascending
if len(results) >= 3:
    # Fit: c1(delta) = c1_0 + c2 * delta^2  (quadratic in delta, symmetric expansion)
    deltas_arr = np.array([r[0] for r in results])
    c1_arr = np.array([r[1] for r in results])

    # Polynomial fit in delta^2
    x = deltas_arr**2
    # Use smallest-delta data points only (to minimize higher-order contamination)
    idx = np.argsort(deltas_arr)[:6]  # smallest 6
    coef = np.polyfit(x[idx], c1_arr[idx], 1)
    c1_extrapolated = coef[1]       # y-intercept (delta=0)
    slope = coef[0]
    print(f"\n  Linear fit c1(delta^2) = c1_0 + c2*delta^2:")
    print(f"    c1_0 (extrapolated) = {c1_extrapolated:.10f}")
    print(f"    c2 slope           = {slope:+.6e}")

    # Compare with hypothesis
    diff = c1_extrapolated - C1_HYPOTHESIS
    rel = abs(diff) / C1_HYPOTHESIS
    print(f"\n  Porownanie z hipoteza 1 - ln(3)/4:")
    print(f"    c1_extrap   = {c1_extrapolated:.10f}")
    print(f"    1 - ln(3)/4 = {C1_HYPOTHESIS:.10f}")
    print(f"    diff        = {diff:+.4e}")
    print(f"    rel err     = {rel:.2e}")

    # Test
    if rel < 1e-4:
        print(f"\n  *** STRONG HIT: c1 = 1 - ln(3)/4  (rel err {rel:.1e}) ***")
        print(f"  *** Implikacja: exp(4(1-c1)) = 3 <=> ln(3) pochodzi z N=3 generacji ***")
    elif rel < 1e-3:
        print(f"\n  POSSIBLE HIT: c1 ≈ 1 - ln(3)/4  (rel err {rel:.1e})")
        print(f"  (Blisko granicy precyzji pomiaru, potrzeba wyzszej rozdzielczosci)")
    else:
        print(f"\n  REJECTED: c1 != 1 - ln(3)/4  (rel err {rel:.1e})")

    # Second hypothesis: exp(4*(1-c1)) = 3 exactly
    exp_check = math.exp(4.0 * (1.0 - c1_extrapolated))
    print(f"\n  Test rownowazny: exp(4*(1-c1_extrap)) = {exp_check:.10f}  (target: 3)")
    print(f"    diff = {exp_check - 3:+.4e}")

    # Try other candidates
    print(f"\n  Sprawdzenie innych kandydatow zamknietych:")
    other_cands = {
        "1 - ln(3)/4":          1 - math.log(3)/4,
        "2/3 + 1/17":           2/3 + 1/17,
        "1/sqrt(2) + 1/50":     1/math.sqrt(2) + 1/50,
        "3/4 - ln(3)/40":       3/4 - math.log(3)/40,
        "(1 + ln(2))/e":        (1 + math.log(2))/math.e,
        "phi/sqrt(5)":          1.6180339887/math.sqrt(5),
    }
    for name, val in sorted(other_cands.items(), key=lambda x: abs(x[1]-c1_extrapolated)):
        d = val - c1_extrapolated
        print(f"    {name:25s} = {val:.10f}  diff = {d:+.4e}")

else:
    print("Za malo punktow do extrapolacji.")
