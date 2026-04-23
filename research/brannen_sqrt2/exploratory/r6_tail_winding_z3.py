#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_tail_winding_z3.py -- Sciezka E: Czy liczba nawiniec n(g0) daje Z3 constraint?

KONTEKST:
  Po negatywnym wyniku Z3 na fazach tail (r6_tail_phase_z3.py) i negatywnym
  wyniku z zasady wariacyjnej (r6_koide_variational.py), pozostaje hipoteza
  TOPOLOGICZNA: moze liczba NAWINIEC (winding number) trojki leptonow
  spelnia n_e + n_mu + n_tau = 0 (mod 3)?

DEFINICJA (winding number w plaszczyznie fazowej):
  theta(r) = atan2(g'(r), g(r) - 1)
  n(g0) = (theta(r_max) - theta(0)) / (2pi)

  To LICZBA CALKOWITA (po unwrap), liczaca ile razy trajektoria (g-1, g')
  okrazyla punkt stabilnej prozni (1, 0) przed osiagnieciem tailu asymptotycznego.

HIPOTEZA E (do przetestowania):
  n_e + n_mu + n_tau = 0 (mod 3) dla fizycznych g0^e=0.869, g0^mu=1.407, g0^tau=1.729

  Jesli TAK -> Z3 dziala na winding numbers -> Koide K=2/3 ma topologiczna derywacje.
  Jesli NIE -> Sciezka E wykluczona.

DODATKOWE TESTY:
  - Skan n(g0) w fizycznym zakresie g0 in (0.5, 2.2)
  - Gdzie sa skoki dyskretne n -> n+1?
  - Czy g0^tau = 1.729 jest SZCZEGOLNYM punktem (granica miedzy skokami)?
  - Czy Hipoteza E jest stabilna wzgledem malych zaburzen g0?

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
    print(f"  {mark}  {name:42s}  {detail}")

# Physical g0 values (from TGP research)
G0_E   = 0.86941
G0_MU  = PHI * G0_E  # phi-drabinka
G0_TAU = 1.7293      # z Koide K=2/3

# PDG masses
M_E   = 0.510999
M_MU  = 105.6584
M_TAU = 1776.86

def rhs(r, y):
    g, gp = y
    if g < 1e-10:
        g = 1e-10
    if r < 1e-12:
        gpp = (1 - g) / 4.0
    else:
        gpp = (1 - g) - (1.0/g)*gp**2 - (2.0/r)*gp
    return [gp, gpp]

def solve_substrate(g0, r_max=400.0, n_points=80000):
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='DOP853', t_eval=r_eval,
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol

def winding_number(r, g, gp, r_cutoff=None):
    """
    Oblicz liczbe nawiniec trajektorii (g-1, g') wokol (0,0).

    theta(r) = atan2(gp, g-1)   (kat biegunowy)
    Po unwrap, n = (theta(r_end) - theta(r_start)) / (2pi)

    r_cutoff: jesli podane, obliczamy tylko do tego punktu (caly soliton zanim
              tail asymptotyczny zdominuje). None = caly przedzial.
    """
    if r_cutoff is not None:
        mask = r <= r_cutoff
        r_f, g_f, gp_f = r[mask], g[mask], gp[mask]
    else:
        r_f, g_f, gp_f = r, g, gp

    # theta(r) = phase angle in (g-1, gp) plane
    theta = np.arctan2(gp_f, g_f - 1.0)
    # Unwrap to get continuous phase
    theta_unwrap = np.unwrap(theta)
    # Total winding (signed, positive = CCW, negative = CW)
    dtheta = theta_unwrap[-1] - theta_unwrap[0]
    n = dtheta / (2.0 * math.pi)
    return n, theta_unwrap

def count_zero_crossings(r, g, r_min=0.1, r_max=400.0):
    """
    Zlicz liczbe przejsc g(r) przez 1 (zero-crossings g-1) w przedziale.
    To alternatywna miara "topologiczna" solitonu.
    """
    mask = (r >= r_min) & (r <= r_max)
    g_f = g[mask]
    # Zero crossings
    sign = np.sign(g_f - 1.0)
    crossings = np.sum(np.abs(np.diff(sign)) > 0)
    return int(crossings)

# ================================================================
print("=" * 72)
print("  R6 Sciezka E: WINDING NUMBER n(g0) -- Z3 constraint na topologii?")
print("=" * 72)

print(f"\n  PHI     = {PHI:.6f}")
print(f"  g0^e    = {G0_E:.5f}")
print(f"  g0^mu   = {G0_MU:.5f} = phi * g0^e")
print(f"  g0^tau  = {G0_TAU:.5f} (z Koide K=2/3)")

# ================================================================
# SEKCJA 1: Oblicz winding n(g0) dla trzech leptonow
# ================================================================
print(f"\n{'='*72}")
print("  1. WINDING NUMBERS DLA e, mu, tau")
print("="*72)

# Dla stabilnosci liczymy winding do r_cutoff gdzie tail asymptotyczny
# juz dominuje ale amplituda |g-1| jeszcze czytelna (np. r = 200).
R_CUTOFF = 200.0

windings = {}
zero_crossings = {}
for name, g0 in [("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)]:
    sol = solve_substrate(g0, r_max=400.0)
    if not sol.success:
        print(f"  {name}: FAILED to solve")
        continue
    n, theta = winding_number(sol.t, sol.y[0], sol.y[1], r_cutoff=R_CUTOFF)
    zc = count_zero_crossings(sol.t, sol.y[0], r_min=0.1, r_max=R_CUTOFF)
    windings[name] = n
    zero_crossings[name] = zc

print(f"\n  {'Lepton':>8s}  {'g0':>10s}  {'n_winding':>14s}  {'n_rounded':>10s}  {'#zero_cross':>12s}")
print(f"  {'-'*8}  {'-'*10}  {'-'*14}  {'-'*10}  {'-'*12}")
for name in ["e", "mu", "tau"]:
    g0 = {"e": G0_E, "mu": G0_MU, "tau": G0_TAU}[name]
    n = windings[name]
    zc = zero_crossings[name]
    print(f"  {name:>8s}  {g0:10.5f}  {n:14.6f}  {round(n):10d}  {zc:12d}")

# ================================================================
# SEKCJA 2: HIPOTEZA E: n_e + n_mu + n_tau = 0 (mod 3)?
# ================================================================
print(f"\n{'='*72}")
print("  2. HIPOTEZA E: Czy n_e + n_mu + n_tau = 0 (mod 3)?")
print("="*72)

# Test dla wartosci zaokraglonych (topologiczne)
n_e = round(windings["e"])
n_mu = round(windings["mu"])
n_tau = round(windings["tau"])
sum_n = n_e + n_mu + n_tau
mod3 = sum_n % 3

print(f"\n  n_e + n_mu + n_tau = {n_e} + {n_mu} + {n_tau} = {sum_n}")
print(f"  Sum mod 3           = {mod3}")

check("Sum of windings equiv 0 (mod 3)", mod3 == 0,
      f"sum = {sum_n}, mod 3 = {mod3}")

# Test dla zero-crossings (alternatywna miara)
zc_e = zero_crossings["e"]
zc_mu = zero_crossings["mu"]
zc_tau = zero_crossings["tau"]
sum_zc = zc_e + zc_mu + zc_tau
mod3_zc = sum_zc % 3

print(f"\n  #zero_cross (e)   = {zc_e}")
print(f"  #zero_cross (mu)  = {zc_mu}")
print(f"  #zero_cross (tau) = {zc_tau}")
print(f"  Sum               = {sum_zc}, mod 3 = {mod3_zc}")

check("Zero-crossings sum equiv 0 (mod 3)", mod3_zc == 0,
      f"sum = {sum_zc}, mod 3 = {mod3_zc}")

# ================================================================
# SEKCJA 3: SKAN n(g0) vs g0 -- gdzie sa skoki dyskretne?
# ================================================================
print(f"\n{'='*72}")
print("  3. SKAN n(g0) -- lokalizacja skokow topologicznych")
print("="*72)

g0_scan = np.linspace(0.3, 2.2, 40)
scan_results = []

print(f"\n  {'g0':>8s}  {'n_winding':>12s}  {'n_int':>6s}  {'#zc':>4s}")
print(f"  {'-'*8}  {'-'*12}  {'-'*6}  {'-'*4}")

for g0 in g0_scan:
    sol = solve_substrate(g0, r_max=300.0, n_points=60000)
    if not sol.success:
        print(f"  {g0:8.4f}  {'FAILED':>12s}")
        continue
    n, _ = winding_number(sol.t, sol.y[0], sol.y[1], r_cutoff=R_CUTOFF)
    zc = count_zero_crossings(sol.t, sol.y[0], r_min=0.1, r_max=R_CUTOFF)
    scan_results.append((g0, n, round(n), zc))
    marker = ""
    # Physical markers
    if abs(g0 - G0_E) < 0.03:
        marker = "  <-- g0^e"
    elif abs(g0 - G0_MU) < 0.03:
        marker = "  <-- g0^mu"
    elif abs(g0 - G0_TAU) < 0.03:
        marker = "  <-- g0^tau"
    print(f"  {g0:8.4f}  {n:12.6f}  {round(n):6d}  {zc:4d}{marker}")

# Znajdz skoki
print(f"\n  SKOKI topologiczne n(g0):")
prev_n = None
for g0, n, ni, zc in scan_results:
    if prev_n is not None and ni != prev_n:
        print(f"    g0 ≈ {g0:.4f}: n {prev_n} -> {ni}")
    prev_n = ni

# ================================================================
# SEKCJA 4: Stabilnosc n wokol g0^tau
# ================================================================
print(f"\n{'='*72}")
print("  4. STABILNOSC n(g0^tau) -- czy 1.729 jest SZCZEGOLNYM punktem?")
print("="*72)

g0_tau_scan = np.linspace(G0_TAU - 0.05, G0_TAU + 0.05, 11)
print(f"\n  {'g0':>8s}  {'n_winding':>12s}  {'n_int':>6s}")
print(f"  {'-'*8}  {'-'*12}  {'-'*6}")
n_tau_values = []
for g0 in g0_tau_scan:
    sol = solve_substrate(g0, r_max=300.0, n_points=60000)
    if not sol.success:
        continue
    n, _ = winding_number(sol.t, sol.y[0], sol.y[1], r_cutoff=R_CUTOFF)
    marker = "  <-- g0^tau" if abs(g0 - G0_TAU) < 0.005 else ""
    print(f"  {g0:8.4f}  {n:12.6f}  {round(n):6d}{marker}")
    n_tau_values.append(round(n))

n_tau_stable = (len(set(n_tau_values)) == 1)
check("n(g0) stable in +/-0.05 around g0^tau", n_tau_stable,
      f"values: {sorted(set(n_tau_values))}")

# ================================================================
# SEKCJA 5: PODSUMOWANIE
# ================================================================
print(f"\n{'='*72}")
print(f"  PODSUMOWANIE: {PASS} PASS, {FAIL} FAIL z {PASS+FAIL} testow")
print("="*72)

if mod3 == 0:
    print("\n  POZYTYWNY WYNIK: Z3 constraint na winding numbers POTWIERDZONY!")
    print("    Nastepny krok: udowodnic topologicznie ze constraint wymusza K=2/3.")
else:
    print(f"\n  NEGATYWNY WYNIK: n_e+n_mu+n_tau = {sum_n} (mod 3 = {mod3}), ")
    print(f"    Z3 constraint nie zachodzi -> Sciezka E wykluczona dla tej definicji.")
    print(f"    Moze nalezy zdefiniowac winding inaczej (np. tail-wrap przed cutoff)?")

print("\n  Windings (round): n_e = {}, n_mu = {}, n_tau = {}".format(n_e, n_mu, n_tau))
print("  Zero crossings:   e = {}, mu = {}, tau = {}".format(zc_e, zc_mu, zc_tau))
