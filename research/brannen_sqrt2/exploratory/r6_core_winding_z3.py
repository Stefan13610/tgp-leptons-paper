#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_core_winding_z3.py -- Sciezka E (v2): CORE winding number

PROBLEM Z v1 (r6_tail_winding_z3.py):
  Winding na calym przedziale [0, R] daje ten sam n = -32 dla wszystkich
  leptonow, bo asymptotyczny tail dominuje (sin(r+delta)/r - ta sama faza r).
  Z3 (mod 3) spelnione trywialnie, nie informatywnie.

NOWA DEFINICJA: CORE WINDING
  Wyodrebnij asymptotyczny tail g_tail(r) = 1 + A*sin(r+delta)/r
  i policzmy winding TYLKO na rezyduum g_core(r) = g(r) - g_tail(r).

  Albo prostsza alternatywa: "lokalny" winding = ile razy (g-1, g') okrazyla
  (0,0) ZANIM wejdzie do asymptotycznego regimu.

ALTERNATYWA: NODE COUNTING na zrenormalizowanym profilu
  Zamiast zliczac wszystkie zera g-1 (ktore rosna jak ~r/pi w tail),
  zliczmy zera g-1 pomniejszonego przez obwiednie OGONA:
  h(r) = (g(r) - 1) * r - A*sin(r + delta)
  Zera h = miejsca gdzie profil ODCHYLA sie od obwiedni asymptotycznej.
  Core-nodes = liczba takich odchylen zanim rezyduum zaniknie.

HIPOTEZA: n_core(g0^e) + n_core(g0^mu) + n_core(g0^tau) = 0 (mod 3)

Author: Claudian
Date: 2026-04-16
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
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

def extract_tail_params(r, g, r_min=100.0, r_max=250.0):
    """Asymptotyczne: u(r) = (g-1)*r = A*sin(r+delta) = c1*sin(r) + c2*cos(r)"""
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.sin(r_f), np.cos(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, u_f, rcond=None)
    c1, c2 = coef
    A = math.sqrt(c1*c1 + c2*c2)
    delta = math.atan2(c2, c1)
    return A, delta

def core_residual(r, g, A_tail, delta_tail):
    """
    Rezyduum po odjeciu obwiedni asymptotycznej:
    h(r) = (g(r) - 1) - A*sin(r+delta)/r    [wszystko w obszarze r > 0]
    Powinno zanikac szybciej niz 1/r w regionie core->tail.
    """
    h = np.zeros_like(g)
    mask = r > 1e-6
    h[mask] = (g[mask] - 1.0) - A_tail * np.sin(r[mask] + delta_tail) / r[mask]
    return h

def count_core_zeros(r, h, r_min=0.5, amplitude_threshold=0.01):
    """
    Policz zera rezyduum h(r) w obszarze gdzie |h| >= threshold.
    Za punkt odciecia uznajemy r_end = maksymalne r dla ktorego
    ruchoma srednia |h| spada ponizej threshold * max|h|.
    """
    mask = r >= r_min
    r_f = r[mask]
    h_f = h[mask]

    # Envelope: zajmuje wartosci po oknie
    abs_h = np.abs(h_f)
    window = 100
    envelope = np.zeros_like(abs_h)
    for i in range(len(abs_h)):
        lo = max(0, i - window)
        hi = min(len(abs_h), i + window)
        envelope[i] = np.max(abs_h[lo:hi])

    # Znajdz gdzie envelope spada ponizej threshold * max
    max_env = np.max(envelope)
    cutoff_idx = len(envelope) - 1
    for i in range(len(envelope)-1, 0, -1):
        if envelope[i] > amplitude_threshold * max_env:
            cutoff_idx = i
            break

    # Liczenie zero-crossings w obszarze core (do cutoff)
    h_core = h_f[:cutoff_idx]
    sign = np.sign(h_core)
    sign[sign == 0] = 1
    crossings = np.sum(np.abs(np.diff(sign)) > 0)
    r_end = r_f[cutoff_idx]
    return int(crossings), r_end

def core_winding(r, g, gp, A_tail, delta_tail, r_min=0.1, r_max_effective=50.0):
    """
    Oblicz winding tylko w regionie CORE (r < r_max_effective gdzie
    profil znaczaco odchyla sie od asymptotycznego tailu).

    r_max_effective dobierane tak aby byc PRZED zdominowaniem przez tail.
    """
    # Residuum w plaszczyznie fazowej: (g - g_tail, g' - g_tail')
    # g_tail = 1 + A*sin(r+delta)/r
    # g_tail' = A*(cos(r+delta)/r - sin(r+delta)/r^2)
    mask = (r >= r_min) & (r <= r_max_effective)
    r_f = r[mask]
    g_f = g[mask]
    gp_f = gp[mask]
    g_tail = 1.0 + A_tail * np.sin(r_f + delta_tail) / r_f
    gp_tail = A_tail * (np.cos(r_f + delta_tail) / r_f - np.sin(r_f + delta_tail) / r_f**2)

    dg = g_f - g_tail
    dgp = gp_f - gp_tail

    theta = np.arctan2(dgp, dg)
    theta_uw = np.unwrap(theta)
    n = (theta_uw[-1] - theta_uw[0]) / (2.0 * math.pi)
    return n

print("=" * 76)
print("  R6 Sciezka E v2: CORE WINDING n(g0) -- Z3 constraint (poprawiona wersja)")
print("=" * 76)

# ================================================================
# SEKCJA 1: Parametry tailu i core residuum dla e, mu, tau
# ================================================================
print(f"\n{'='*76}")
print("  1. EKSTRAKCJA (A_tail, delta_tail) + core residuum")
print("="*76)

results = {}
print(f"\n  {'Lepton':>8s}  {'g0':>9s}  {'A_tail':>12s}  {'delta':>10s}  {'#core_zeros':>12s}  {'r_end':>7s}")
print(f"  {'-'*8}  {'-'*9}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*7}")

for name, g0 in [("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)]:
    sol = solve_substrate(g0, r_max=400.0)
    A, delta = extract_tail_params(sol.t, sol.y[0])
    h = core_residual(sol.t, sol.y[0], A, delta)
    zeros_core, r_end = count_core_zeros(sol.t, h, r_min=0.5, amplitude_threshold=0.02)
    results[name] = {
        "g0": g0, "A": A, "delta": delta,
        "core_zeros": zeros_core, "r_end": r_end,
        "sol": sol, "h": h
    }
    print(f"  {name:>8s}  {g0:9.5f}  {A:12.8f}  {delta:10.5f}  {zeros_core:12d}  {r_end:7.2f}")

# ================================================================
# SEKCJA 2: HIPOTEZA E: n_core(e) + n_core(mu) + n_core(tau) = 0 (mod 3)?
# ================================================================
print(f"\n{'='*76}")
print("  2. HIPOTEZA E: Z3 na core_zeros ?")
print("="*76)

n_e = results["e"]["core_zeros"]
n_mu = results["mu"]["core_zeros"]
n_tau = results["tau"]["core_zeros"]
sum_n = n_e + n_mu + n_tau
mod3 = sum_n % 3

print(f"\n  n_core(e)   = {n_e}")
print(f"  n_core(mu)  = {n_mu}")
print(f"  n_core(tau) = {n_tau}")
print(f"  SUM         = {sum_n}")
print(f"  SUM mod 3   = {mod3}")

check("Nontrivial: not all equal", len({n_e, n_mu, n_tau}) > 1,
      f"values {{e={n_e}, mu={n_mu}, tau={n_tau}}}")
check("SUM equiv 0 (mod 3)", mod3 == 0,
      f"sum={sum_n}, mod 3={mod3}")

# ================================================================
# SEKCJA 3: Core winding w plaszczyznie fazowej (rezyduum)
# ================================================================
print(f"\n{'='*76}")
print("  3. CORE WINDING w plaszczyznie (dg, dgp) = residuum")
print("="*76)

print(f"\n  {'Lepton':>8s}  {'g0':>9s}  {'core_winding':>16s}  {'rounded':>8s}")
print(f"  {'-'*8}  {'-'*9}  {'-'*16}  {'-'*8}")

windings = {}
for name in ["e", "mu", "tau"]:
    r = results[name]["sol"].t
    g = results[name]["sol"].y[0]
    gp = results[name]["sol"].y[1]
    A = results[name]["A"]
    delta = results[name]["delta"]
    # Wybierz r_max_effective w oparciu o r_end
    r_max_eff = min(results[name]["r_end"] + 5.0, 50.0)
    nw = core_winding(r, g, gp, A, delta,
                      r_min=0.5, r_max_effective=r_max_eff)
    windings[name] = nw
    print(f"  {name:>8s}  {results[name]['g0']:9.5f}  {nw:16.6f}  {round(nw):8d}")

nw_e = round(windings["e"])
nw_mu = round(windings["mu"])
nw_tau = round(windings["tau"])
sum_nw = nw_e + nw_mu + nw_tau
mod3_nw = sum_nw % 3
print(f"\n  SUM = {sum_nw}, mod 3 = {mod3_nw}")
check("Core winding SUM equiv 0 (mod 3)", mod3_nw == 0,
      f"sum={sum_nw}, mod 3={mod3_nw}")
check("Core winding NONTRIVIAL (not all equal)", len({nw_e,nw_mu,nw_tau}) > 1,
      f"{{e={nw_e}, mu={nw_mu}, tau={nw_tau}}}")

# ================================================================
# SEKCJA 4: Skan core_zeros vs g0 w szerszym zakresie
# ================================================================
print(f"\n{'='*76}")
print("  4. SKAN n_core(g0) -- identyfikacja plateau")
print("="*76)

g0_scan = np.linspace(0.5, 2.2, 35)
print(f"\n  {'g0':>8s}  {'A_tail':>10s}  {'delta':>8s}  {'n_core':>7s}  {'r_end':>7s}")
print(f"  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*7}  {'-'*7}")

for g0 in g0_scan:
    sol = solve_substrate(g0, r_max=400.0, n_points=60000)
    A, d = extract_tail_params(sol.t, sol.y[0])
    h = core_residual(sol.t, sol.y[0], A, d)
    n_c, r_end = count_core_zeros(sol.t, h, r_min=0.5, amplitude_threshold=0.02)
    marker = ""
    if abs(g0 - G0_E) < 0.03:
        marker = "  <-- g0^e"
    elif abs(g0 - G0_MU) < 0.03:
        marker = "  <-- g0^mu"
    elif abs(g0 - G0_TAU) < 0.03:
        marker = "  <-- g0^tau"
    print(f"  {g0:8.4f}  {A:10.6f}  {d:8.4f}  {n_c:7d}  {r_end:7.2f}{marker}")

# ================================================================
print(f"\n{'='*76}")
print(f"  PODSUMOWANIE: {PASS} PASS, {FAIL} FAIL z {PASS+FAIL} testow")
print("="*76)

if PASS == PASS + FAIL:
    print("\n  MOZLIWY POZYTYWNY WYNIK -- dalsza weryfikacja wymagana")
else:
    print(f"\n  MIESZANY WYNIK ({PASS}/{PASS+FAIL}) -- interpretacja:")
    if not (len({n_e, n_mu, n_tau}) > 1):
        print("    Brak zroznicowania miedzy leptonami - trzeba redefinicji.")
    if mod3 != 0:
        print("    Z3 constraint wykluczone dla core_zeros.")
