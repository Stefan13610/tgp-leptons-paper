#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c1_richardson.py
===================

Polepszona extrapolacja c1 z dodatkowymi rzedami.
Sprawdzam czy O(delta^4) zmienia c1_0.
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import math

# Dane z r6_c1_high_precision.py (skopiowane rekordy)
data = [
    (0.05000, 0.7251350584),
    (0.03000, 0.7252137877),
    (0.02000, 0.7252383815),
    (0.01500, 0.7252469895),
    (0.01000, 0.7252531479),
    (0.00800, 0.7252549049),
    (0.00600, 0.7252562739),
    (0.00400, 0.7252572234),
    (0.00200, 0.7252577354),
]

deltas = np.array([d for d, _ in data])
c1_arr = np.array([c for _, c in data])

print("=" * 72)
print("  Richardson extrapolacja c1(delta) -> c1(0)")
print("=" * 72)
print()
print(f"  Surowe dane (9 probek, delta = 0.002 do 0.05):")
for d, c in data:
    print(f"    delta = {d:.5f}   c1 = {c:.10f}")

# Try different polynomial fits
print(f"\n{'-'*72}")
print("  Fit polinomial c1(delta) = a0 + a2*delta^2 + a4*delta^4")
print(f"{'-'*72}")

# Fit 1: quadratic (just delta^2)
A = np.column_stack([np.ones_like(deltas), deltas**2])
coef, res, *_ = np.linalg.lstsq(A, c1_arr, rcond=None)
print(f"\n  Linear in delta^2:   c1_0 = {coef[0]:.10f}  + ({coef[1]:+.4e})*delta^2")
resid = c1_arr - A @ coef
print(f"    max residual = {np.max(np.abs(resid)):.2e}")

# Fit 2: quadratic + quartic
A2 = np.column_stack([np.ones_like(deltas), deltas**2, deltas**4])
coef2, *_ = np.linalg.lstsq(A2, c1_arr, rcond=None)
print(f"\n  With delta^4:        c1_0 = {coef2[0]:.10f}  + ({coef2[1]:+.4e})*delta^2  + ({coef2[2]:+.4e})*delta^4")
resid2 = c1_arr - A2 @ coef2
print(f"    max residual = {np.max(np.abs(resid2)):.2e}")

# Fit 3: use only smallest deltas (more accurate for extrapolation)
small_idx = deltas <= 0.015
A3 = np.column_stack([np.ones(small_idx.sum()), deltas[small_idx]**2])
coef3, *_ = np.linalg.lstsq(A3, c1_arr[small_idx], rcond=None)
print(f"\n  Only smallest 5 pts: c1_0 = {coef3[0]:.10f}  + ({coef3[1]:+.4e})*delta^2")

# Compare with candidates
print(f"\n{'-'*72}")
print("  Porownanie z kandydatami zamknietymi")
print(f"{'-'*72}")

c1_best = coef2[0]  # use best estimate (with delta^4 correction)
print(f"\n  Najlepsza ekstrapolacja: c1_0 = {c1_best:.10f}")

candidates = {
    "1 - ln(3)/4":                1 - math.log(3)/4,
    "1 - ln(3)/4 - 1e-4":         1 - math.log(3)/4 - 1.0e-4,
    "(1 + 1/25)*0.7":             (1 + 1/25)*0.7,   # 0.728 - probably off
    "2/3 + 1/(4*pi)":             2/3 + 1/(4*math.pi),
    "2/3 + ln(e^pi-20)/100":      2/3 + math.log(math.exp(math.pi)-20)/100,  # empirical
    "sqrt(pi)/sqrt(2*pi - 1)":    math.sqrt(math.pi)/math.sqrt(2*math.pi - 1),
    "(pi-sqrt(2))/(pi-1/2)":      (math.pi-math.sqrt(2))/(math.pi-0.5),
    "1 - 1/sqrt(2*pi)":           1 - 1/math.sqrt(2*math.pi),
    "2 - 2/sqrt(pi^2+1)":         2 - 2/math.sqrt(math.pi**2+1),
    "pi/(2*phi)":                 math.pi/(2*1.6180339887),
    "1 - 1/(pi*sqrt(2))":         1 - 1/(math.pi*math.sqrt(2)),
    "3/4 * (1 - delta0)":         None,  # placeholder
}
del candidates["3/4 * (1 - delta0)"]

print(f"\n  {'Candidate':35s} {'Value':>14s}  {'diff':>12s}  {'rel err':>10s}")
print(f"  {'-'*35} {'-'*14}  {'-'*12}  {'-'*10}")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-c1_best)):
    diff = val - c1_best
    rel = abs(diff)/c1_best
    print(f"  {name:35s} {val:14.10f}  {diff:+12.4e}  {rel:10.2e}")

# Final conclusion
print(f"\n{'-'*72}")
print("  Wnioski")
print(f"{'-'*72}")
tol_hypothesis = 1e-5
hyp = 1 - math.log(3)/4
err_hyp = abs(hyp - c1_best) / c1_best
if err_hyp < tol_hypothesis:
    print(f"\n  *** c1 = 1 - ln(3)/4 POTWIERDZONE do {tol_hypothesis:.0e} ***")
else:
    print(f"\n  c1 != 1 - ln(3)/4  na poziomie 1e-5")
    print(f"  (residual {err_hyp:.2e}, tolerancja {tol_hypothesis:.0e})")
    print(f"\n  Wyznaczona wartosc c1 = {c1_best:.8f}")
    print(f"  Brak oczywistego zamknietego kandydatu z listy")
    print(f"  Potencjalnie c1 zalezne od c=1 lub (1 - 1/r_typical) albo")
    print(f"  subtelnego numeryka rozwiazania ODE")

# Dodatkowe sprawdzenia: c1 * (stale)
print(f"\n  Przeciwne: rozne kombinacje z c1:")
print(f"    1/c1          = {1/c1_best:.8f}")
print(f"    c1 * pi       = {c1_best * math.pi:.8f}")
print(f"    c1 * e        = {c1_best * math.e:.8f}")
print(f"    c1 * sqrt(2)  = {c1_best * math.sqrt(2):.8f}")
print(f"    c1 / 0.5      = {c1_best / 0.5:.8f}")
print(f"    1 - c1        = {1 - c1_best:.8f}")
print(f"    4 * (1-c1)    = {4*(1-c1_best):.8f}   (ln 3 = {math.log(3):.8f})")
