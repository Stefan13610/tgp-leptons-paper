#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_extrema_index_z3.py -- Sciezka E (v3): TURNING NUMBER (liczba lokalnych ekstremow)

Po negatywnych wynikach:
  v1: total winding = -32 dla wszystkich (tail dominuje)
  v2: core_zeros = 0 dla wszystkich (wszystkie groundstate solitonami)

Inna wielkosc topologiczna, ktora MOZE roznicowac leptony:

TURNING NUMBER T(g0) = liczba zmian znaku g'(r) w CALYM przedziale [0, R_max].
Czyli liczba lokalnych ekstremow g(r).

Dla soliton bez wezlow wewnetrznych, T rosnie zgodnie z liczba oscylacji
ogona. Roznicuje sie miedzy leptonami przez dlugosc wyjscia z core
(im glebszy soliton, tym pozniej wchodzi w oscylujacy tail).

ALTERNATYWA: BOUNCE NUMBER = liczba razy g'(r) zmienia znak w CORE
  (przed pierwszym przejsciem przez 1).
  Dla e/mu/tau: g(0)=g0 rozne, g'(0)=0, g'(r) zmienia znak gdy g(r) osiagnie
  dolny zwrot trajektorii fazowej. Mozliwe zroznicowanie.

HIPOTEZA (weak): T(g0^e) + T(g0^mu) + T(g0^tau) = 0 (mod 3)

Author: Claudian
Date: 2026-04-16
"""

import numpy as np
from scipy.integrate import solve_ivp
import math

PHI = (1 + math.sqrt(5)) / 2

PASS = 0
FAIL = 0

def check(name, cond, detail=""):
    global PASS, FAIL
    mark = "PASS" if cond else "FAIL"
    if cond:
        PASS += 1
    else:
        FAIL += 1
    print(f"  {mark}  {name:46s}  {detail}")

G0_E   = 0.86941
G0_MU  = PHI * G0_E
G0_TAU = 1.7293

def rhs(r, y):
    g, gp = y
    if g < 1e-10: g = 1e-10
    if r < 1e-12: gpp = (1 - g) / 4.0
    else: gpp = (1 - g) - (1.0/g)*gp**2 - (2.0/r)*gp
    return [gp, gpp]

def solve_substrate(g0, r_max=200.0, n_points=80000):
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol

def count_extrema(r, g, r_max_window):
    """Liczba lokalnych ekstremow g(r) w [0, r_max_window]."""
    mask = r <= r_max_window
    g_f = g[mask]
    dg = np.diff(g_f)
    sign_dg = np.sign(dg)
    # Number of sign changes of derivative
    changes = np.sum(np.abs(np.diff(sign_dg)) > 0)
    return int(changes)

def count_crossings_of_one(r, g, r_max_window):
    """Liczba przejsc g=1 w [0, r_max_window]."""
    mask = r <= r_max_window
    g_f = g[mask]
    sign = np.sign(g_f - 1.0)
    sign[sign == 0] = 1
    return int(np.sum(np.abs(np.diff(sign)) > 0))

def first_crossing_of_one(r, g):
    """r gdzie g po raz pierwszy przechodzi przez 1 (konczy core region)."""
    g0 = g[0]
    target = 1.0
    # Jesli g0 > 1: szukamy kiedy g <= 1
    # Jesli g0 < 1: szukamy kiedy g >= 1
    if g0 > 1:
        mask = g <= target
    else:
        mask = g >= target
    idx = np.argmax(mask)
    if not mask[idx]:
        return None
    return r[idx]

def bounce_count_in_core(r, g, gp, r_max_window):
    """Liczba razy gp=0 w przedziale [0, r_max_window]."""
    mask = r <= r_max_window
    gp_f = gp[mask]
    sign = np.sign(gp_f)
    sign[sign == 0] = 1
    return int(np.sum(np.abs(np.diff(sign)) > 0))

print("=" * 76)
print("  R6 Sciezka E v3: TURNING NUMBER / EXTREMA COUNT -- Z3 test")
print("=" * 76)

# ================================================================
# SEKCJA 1: Oblicz rozne wskazniki topologiczne
# ================================================================
print(f"\n{'='*76}")
print("  1. WSKAZNIKI TOPOLOGICZNE dla e, mu, tau")
print("="*76)

R_MAX = 30.0  # okienko core (po ktorym tail dominuje)

results = {}
for name, g0 in [("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)]:
    sol = solve_substrate(g0, r_max=200.0)
    r, g, gp = sol.t, sol.y[0], sol.y[1]

    # r_match: r gdzie profil "osiaga" tail (bierzemy pierwsze g=1 + offset)
    r_first = first_crossing_of_one(r, g)
    r_match = r_first if r_first is not None else 10.0

    # Rozne wskazniki
    n_extr = count_extrema(r, g, R_MAX)
    n_cross = count_crossings_of_one(r, g, R_MAX)
    n_bounce = bounce_count_in_core(r, g, gp, R_MAX)
    # Core extrema: tylko do r_first + 2
    r_core = min(r_match + 2.0 if r_first else 8.0, 15.0)
    n_extr_core = count_extrema(r, g, r_core)
    n_bounce_core = bounce_count_in_core(r, g, gp, r_core)

    results[name] = {
        "g0": g0, "r_first": r_first,
        "n_extr_30": n_extr,
        "n_cross_30": n_cross,
        "n_bounce_30": n_bounce,
        "n_extr_core": n_extr_core,
        "n_bounce_core": n_bounce_core,
        "r_core": r_core
    }

print(f"\n  {'Lepton':>8s}  {'g0':>9s}  {'r_first':>8s}  {'#extr(30)':>10s}  {'#cross(30)':>11s}  {'#bnc_core':>10s}  {'#ext_core':>10s}")
print(f"  {'-'*8}  {'-'*9}  {'-'*8}  {'-'*10}  {'-'*11}  {'-'*10}  {'-'*10}")
for name in ["e", "mu", "tau"]:
    r = results[name]
    rf = f"{r['r_first']:.3f}" if r['r_first'] else "---"
    print(f"  {name:>8s}  {r['g0']:9.5f}  {rf:>8s}  {r['n_extr_30']:10d}  {r['n_cross_30']:11d}  {r['n_bounce_core']:10d}  {r['n_extr_core']:10d}")

# ================================================================
# SEKCJA 2: HIPOTEZY Z3
# ================================================================
print(f"\n{'='*76}")
print("  2. HIPOTEZY Z3 -- czy ktoras wielkosc spelnia sum mod 3 = 0?")
print("="*76)

metrics = {
    "n_extr_30": [results[n]["n_extr_30"] for n in ["e","mu","tau"]],
    "n_cross_30": [results[n]["n_cross_30"] for n in ["e","mu","tau"]],
    "n_bounce_core": [results[n]["n_bounce_core"] for n in ["e","mu","tau"]],
    "n_extr_core": [results[n]["n_extr_core"] for n in ["e","mu","tau"]],
}

for metric_name, vals in metrics.items():
    s = sum(vals)
    m = s % 3
    distinct = len(set(vals)) > 1
    print(f"\n  {metric_name}: e={vals[0]}, mu={vals[1]}, tau={vals[2]}, SUM={s}, mod3={m}")
    check(f"{metric_name}: nontrivial differentiation", distinct, f"vals={vals}")
    check(f"{metric_name}: sum equiv 0 (mod 3)", m == 0, f"sum={s}")

# ================================================================
# SEKCJA 3: SKAN metryk vs g0
# ================================================================
print(f"\n{'='*76}")
print("  3. SKAN: jak wskazniki zmieniaja sie z g0?")
print("="*76)

print(f"\n  {'g0':>8s}  {'r_first':>8s}  {'#ext_core':>10s}  {'#bnc_core':>10s}  {'#extr(30)':>10s}  {'#cross(30)':>11s}")
print(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*11}")

g0_scan = np.linspace(0.4, 2.3, 40)
for g0 in g0_scan:
    sol = solve_substrate(g0, r_max=200.0, n_points=60000)
    r, g, gp = sol.t, sol.y[0], sol.y[1]
    r_first = first_crossing_of_one(r, g)
    r_core = min((r_first if r_first else 8.0) + 2.0, 15.0)
    n_ec = count_extrema(r, g, r_core)
    n_bc = bounce_count_in_core(r, g, gp, r_core)
    n_e30 = count_extrema(r, g, R_MAX)
    n_c30 = count_crossings_of_one(r, g, R_MAX)
    rf = f"{r_first:.3f}" if r_first else "---"
    marker = ""
    if abs(g0 - G0_E) < 0.03: marker = "  e"
    elif abs(g0 - G0_MU) < 0.03: marker = "  mu"
    elif abs(g0 - G0_TAU) < 0.03: marker = "  tau"
    print(f"  {g0:8.4f}  {rf:>8s}  {n_ec:10d}  {n_bc:10d}  {n_e30:10d}  {n_c30:11d}{marker}")

# ================================================================
print(f"\n{'='*76}")
print(f"  PODSUMOWANIE: {PASS} PASS, {FAIL} FAIL z {PASS+FAIL} testow")
print("="*76)

if FAIL == 0:
    print("\n  WSZYSTKIE TESTY PASS")
elif PASS >= FAIL:
    print(f"\n  MIESZANE: {PASS} PASS, {FAIL} FAIL")
else:
    print(f"\n  NEGATYWNE: {PASS} PASS, {FAIL} FAIL")
    print("  Wniosek: Z3 constraint na topologicznych indeksach e/mu/tau NIE zachodzi.")
    print("  Wszystkie leptony sa w tym samym 'sektorze topologicznym' (groundstate),")
    print("  wiec Koide K=2/3 NIE WYNIKA z topologii liczby wezlow/ekstremow.")
