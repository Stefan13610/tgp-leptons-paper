#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_tail_phase_z3.py -- Czy fazy delta_i ogonow solitonowych sa Z3-symetryczne?

HIPOTEZA: Tail asymptotyczny g(r) = 1 + A*sin(r+delta)/r dla trzech leptonow
ma fazy delta_e, delta_mu, delta_tau rozlozone na 120 stopni (Z3).

Jesli TAK -> Brannen sqrt(m_i) = a(1 + b*cos(theta + 2pi*i/3)) ma NATURALNE
wyjasnienie jako rzut Z3-rownowaznej struktury fazowej na amplitude A.

OCZEKIWANY MECHANIZM:
  Jesli delta_e = delta_0,  delta_mu = delta_0 + 2pi/3,  delta_tau = delta_0 + 4pi/3
  (mod 2pi) i A_i = A(delta_i) dla pewnej funkcji A(delta),
  to amplitudy A_i automatycznie spelniaja strukture 120-stopniowa.

Ale to jest mocna HIPOTEZA -- testujemy numerycznie.

TEST:
  1. Rozwiaz ODE dla g0^e, g0^mu, g0^tau fizycznych
  2. Ekstraktuj (A_i, delta_i) z dopasowania g(r) - 1 = A*sin(r+delta)/r
  3. Oblicz roznice delta_j - delta_i (mod 2pi)
  4. Porownaj z 2pi/3 = 120 stopni

Dodatkowo:
  - Skan delta(g0) po fizycznym zakresie g0
  - Sprawdzic CIAGLA zmiana delta w g0
  - Czy Brannen A(delta) = a + b*cos(delta) pasuje?

Author: Claudian
Date: 2026-04-16
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize, curve_fit
import math

PHI = (1 + math.sqrt(5)) / 2
TWO_PI_3 = 2 * math.pi / 3

PASS = 0
FAIL = 0

def check(name, cond, detail=""):
    global PASS, FAIL
    mark = "PASS" if cond else "FAIL"
    if cond:
        PASS += 1
    else:
        FAIL += 1
    print(f"  {mark}  {name:40s}  {detail}")

# Physical g0 values (from Koide + phi-ladder inversion, see r6_koide_from_ode.py)
G0_E   = 0.86941
G0_MU  = PHI * G0_E  # phi-drabinka
G0_TAU = 1.7293  # z Koide K=2/3

# PDG masses for verification
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
                    method='RK45', t_eval=r_eval,
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol

def fit_tail_amp_phase(r, g, r_min=80.0, r_max=250.0):
    """
    Dopasuj g(r) - 1 = A*sin(r+delta)/r na okienku (r_min, r_max).
    Wyekstraktuj (A, delta) zamiast tylko A.

    Metoda: u(r) = (g-1)*r ~ A*sin(r+delta) = A*(cos(delta)*sin(r) + sin(delta)*cos(r))
          = c1*sin(r) + c2*cos(r)
    z c1 = A*cos(delta), c2 = A*sin(delta)
    => A = sqrt(c1^2 + c2^2), delta = atan2(c2, c1)
    """
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 20:
        return None, None, None
    u_f = (g[mask] - 1.0) * r_f
    # Dwie kolumny: sin(r) i cos(r)
    X = np.column_stack([np.sin(r_f), np.cos(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, u_f, rcond=None)
    c1, c2 = coef
    A = math.sqrt(c1*c1 + c2*c2)
    delta = math.atan2(c2, c1)  # phase of A*sin(r+delta)
    # Reszty (rmse/A < 5% => dobry fit oscylacyjny)
    y_hat = c1*np.sin(r_f) + c2*np.cos(r_f)
    rmse = float(np.sqrt(np.mean((u_f - y_hat)**2)))
    return A, delta, rmse / max(A, 1e-12)

def phase_diff_mod(a, b, period=2*math.pi):
    """Zwroc roznice (a-b) mod period w zakresie [-period/2, period/2]"""
    d = (a - b) % period
    if d > period/2:
        d -= period
    return d

# ================================================================
print("=" * 70)
print("  R6: FAZY OGONA Z3 -- Czy delta_i sa 120 stopni apart?")
print("=" * 70)

print(f"\n  PHI     = {PHI:.6f}")
print(f"  g0^e    = {G0_E:.5f}")
print(f"  g0^mu   = {G0_MU:.5f} = phi * g0^e")
print(f"  g0^tau  = {G0_TAU:.5f} (z Koide K=2/3)")
print(f"  r21 PDG = {M_MU/M_E:.4f}")
print(f"  r31 PDG = {M_TAU/M_E:.4f}")

# ================================================================
# SECTION 1: Ekstraktuj (A, delta) dla trzech leptonow
# ================================================================
print(f"\n{'='*70}")
print("  1. AMPLITUDY I FAZY OGONA DLA e, mu, tau")
print("="*70)

results = {}
for name, g0 in [("e", G0_E), ("mu", G0_MU), ("tau", G0_TAU)]:
    sol = solve_substrate(g0)
    if not sol.success:
        print(f"  {name}: FAILED to solve")
        continue
    A, delta, rel_err = fit_tail_amp_phase(sol.t, sol.y[0])
    if A is None:
        print(f"  {name}: FAILED to fit tail")
        continue
    results[name] = {
        "g0": g0,
        "A": A,
        "delta": delta,
        "delta_deg": math.degrees(delta),
        "rel_err": rel_err,
        "A4": A**4,
    }

print(f"\n  {'Lepton':>8s}  {'g0':>10s}  {'A_tail':>12s}  {'delta (rad)':>14s}  {'delta (deg)':>14s}  {'rmse/A':>10s}")
print(f"  {'-'*8}  {'-'*10}  {'-'*12}  {'-'*14}  {'-'*14}  {'-'*10}")
for name, r in results.items():
    print(f"  {name:>8s}  {r['g0']:10.5f}  {r['A']:12.8f}  {r['delta']:14.6f}  {r['delta_deg']:14.4f}  {r['rel_err']:10.2e}")

# ================================================================
# SECTION 2: Roznice faz -- Z3?
# ================================================================
print(f"\n{'='*70}")
print("  2. ROZNICE FAZ -- test Z3 (120 stopni)")
print("="*70)

if set(results.keys()) >= {"e", "mu", "tau"}:
    d_e   = results["e"]["delta"]
    d_mu  = results["mu"]["delta"]
    d_tau = results["tau"]["delta"]

    diff_mu_e   = phase_diff_mod(d_mu,  d_e)
    diff_tau_mu = phase_diff_mod(d_tau, d_mu)
    diff_tau_e  = phase_diff_mod(d_tau, d_e)

    print(f"\n  delta_mu  - delta_e   = {diff_mu_e:+.6f} rad = {math.degrees(diff_mu_e):+8.3f} deg")
    print(f"  delta_tau - delta_mu  = {diff_tau_mu:+.6f} rad = {math.degrees(diff_tau_mu):+8.3f} deg")
    print(f"  delta_tau - delta_e   = {diff_tau_e:+.6f} rad = {math.degrees(diff_tau_e):+8.3f} deg")

    print(f"\n  Target Z3: 2pi/3 = {TWO_PI_3:.6f} rad = 120.000 deg")
    print(f"  Abs od 120 stopni:")
    for name, d in [("mu-e", diff_mu_e), ("tau-mu", diff_tau_mu), ("tau-e", diff_tau_e)]:
        # Rozszerz do 120 lub -120 (obie strony cykla)
        dist_120 = min(abs(math.degrees(d) - 120), abs(math.degrees(d) + 120), abs(math.degrees(d) - 240), abs(math.degrees(d) + 240))
        print(f"    {name:10s}: diff = {math.degrees(d):+8.3f} deg, odl od +/-120 = {dist_120:.3f} deg")

    # Najprostszy test: czy trzy fazy pasuja do Z3 (120 stopni apart)?
    phases = [d_e, d_mu, d_tau]
    # Znormalizuj do [0, 2pi)
    phases_norm = [p % (2*math.pi) for p in phases]
    # Sortuj
    phases_sorted = sorted(phases_norm)
    # Oblicz roznice cykliczne
    diffs = []
    for i in range(3):
        next_i = (i+1) % 3
        if next_i == 0:
            d = phases_sorted[0] + 2*math.pi - phases_sorted[2]
        else:
            d = phases_sorted[next_i] - phases_sorted[i]
        diffs.append(d)

    print(f"\n  Fazy sorted: {[math.degrees(p) for p in phases_sorted]}")
    print(f"  Cykliczne roznice (deg): {[math.degrees(d) for d in diffs]}")
    print(f"  Max odchylenie od 120 (deg): {max(abs(math.degrees(d) - 120) for d in diffs):.4f}")

    max_dev = max(abs(math.degrees(d) - 120) for d in diffs)
    check("T1: Roznice faz = 120 stopni (tolerancja 5 deg)",
          max_dev < 5.0,
          f"max dev = {max_dev:.3f} deg")
    check("T2: Roznice faz = 120 stopni (tolerancja 15 deg)",
          max_dev < 15.0,
          f"max dev = {max_dev:.3f} deg")

# ================================================================
# SECTION 3: Skan delta(g0) w zakresie fizycznym
# ================================================================
print(f"\n{'='*70}")
print("  3. SKAN delta(g0) -- jak faza zalezy od g0?")
print("="*70)

g0_scan = np.concatenate([
    np.linspace(0.3, 0.98, 18),
    np.linspace(1.02, 2.15, 25),
])

scan_results = []
for g0 in g0_scan:
    sol = solve_substrate(g0)
    if not sol.success:
        continue
    A, delta, rel_err = fit_tail_amp_phase(sol.t, sol.y[0])
    if A is not None and rel_err < 0.1:
        scan_results.append({"g0": g0, "A": A, "delta": delta, "err": rel_err})

print(f"\n  Ilosc skanow: {len(scan_results)}")
print(f"\n  {'g0':>8s}  {'A':>12s}  {'delta':>12s}  {'delta/pi':>10s}  {'A*cos(d)':>12s}  {'A*sin(d)':>12s}")
print(f"  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*12}")
for s in scan_results:
    marker = ""
    if abs(s["g0"] - G0_E)   < 0.01: marker = " <-- e"
    elif abs(s["g0"] - G0_MU)  < 0.02: marker = " <-- mu"
    elif abs(s["g0"] - G0_TAU) < 0.02: marker = " <-- tau"
    print(f"  {s['g0']:8.4f}  {s['A']:12.6f}  {s['delta']:12.6f}  {s['delta']/math.pi:10.5f}  {s['A']*math.cos(s['delta']):12.6f}  {s['A']*math.sin(s['delta']):12.6f}{marker}")

# ================================================================
# SECTION 4: Czy A vs delta pasuje do Brannen a + b*cos(theta)?
# ================================================================
print(f"\n{'='*70}")
print("  4. A(delta) -- czy pasuje do Brannen sqrt(m) = a + b*cos(theta)?")
print("="*70)

if len(scan_results) > 5 and "e" in results and "mu" in results and "tau" in results:
    # Model Brannena: sqrt(m) ~ A^2 (bo m ~ A^4). So to maj to:
    # A(delta) = sqrt(a^2 + 2ab*cos(delta) + b^2*cos^2(delta))? Not exact cos.
    #
    # Ale prosta hipoteza: A = a0 + b0*cos(delta - delta_0)
    # (Brannen parametryzuje SQRT(m), nie m. Ale SQRT(m) = A^2. So:
    # A^2 = a + b*cos(delta - delta_0).)
    #
    # TESTUJEMY obie formy.

    deltas = np.array([s["delta"] for s in scan_results])
    As = np.array([s["A"] for s in scan_results])

    # Fit A = a + b*cos(delta - delta_0)
    def brannen_A(delta, a, b, delta_0):
        return a + b * np.cos(delta - delta_0)

    # Fit A^2 = a + b*cos(delta - delta_0) (Brannen na sqrt(m) = A^2)
    def brannen_A2(delta, a, b, delta_0):
        return a + b * np.cos(delta - delta_0)

    try:
        popt_A, _ = curve_fit(brannen_A, deltas, As, p0=[np.mean(As), np.std(As), 0])
        A_fit = brannen_A(deltas, *popt_A)
        rmse_A = np.sqrt(np.mean((As - A_fit)**2))
        print(f"\n  Fit A = a + b*cos(delta - delta_0):")
        print(f"    a       = {popt_A[0]:.6f}")
        print(f"    b       = {popt_A[1]:.6f}")
        print(f"    delta_0 = {popt_A[2]:.6f} rad = {math.degrees(popt_A[2]):.3f} deg")
        print(f"    B = b/a = {popt_A[1]/popt_A[0]:.6f} (target sqrt(2) = {math.sqrt(2):.6f})")
        print(f"    RMSE    = {rmse_A:.6f} (relative: {rmse_A/np.mean(np.abs(As))*100:.2f}%)")
    except Exception as e:
        print(f"  Fit A: failed ({e})")

    try:
        popt_A2, _ = curve_fit(brannen_A2, deltas, As**2, p0=[np.mean(As**2), np.std(As**2), 0])
        A2_fit = brannen_A2(deltas, *popt_A2)
        rmse_A2 = np.sqrt(np.mean((As**2 - A2_fit)**2))
        print(f"\n  Fit A^2 = a + b*cos(delta - delta_0) (Brannen na sqrt(m)):")
        print(f"    a       = {popt_A2[0]:.6f}")
        print(f"    b       = {popt_A2[1]:.6f}")
        print(f"    delta_0 = {popt_A2[2]:.6f} rad = {math.degrees(popt_A2[2]):.3f} deg")
        print(f"    B = b/a = {popt_A2[1]/popt_A2[0]:.6f} (target sqrt(2) = {math.sqrt(2):.6f})")
        print(f"    RMSE    = {rmse_A2:.6f} (relative: {rmse_A2/np.mean(np.abs(As**2))*100:.2f}%)")

        B_fit = popt_A2[1] / popt_A2[0]
        check("T3: B = b/a bliskie sqrt(2) (toleracja 20%)",
              abs(B_fit - math.sqrt(2)) / math.sqrt(2) < 0.20,
              f"B_fit = {B_fit:.4f}, target = {math.sqrt(2):.4f}")
    except Exception as e:
        print(f"  Fit A^2: failed ({e})")

# ================================================================
# SECTION 5: Alternatywna hipoteza -- delta dzialek okresu 2pi wobec r_max
# ================================================================
print(f"\n{'='*70}")
print("  5. ALTERNATYWNE: czy g0 -> delta z ciagla krzywa?")
print("="*70)

print("""
  Mozliwosc: delta(g0) jest CIAGLA, a fazy leptonowe NIE sa Z3-rownowazne.
  Zamiast tego mozliwa inna struktura, np. delta zawsze w [-pi, pi]
  i trzy leptony po prostu wybieraja rozne g0 dajace delta na skali liniowej.

  Sprawdzmy: czy delta(g0) jest monotoniczna lub ma specjalna strukture?
""")

if len(scan_results) > 5:
    deltas_scan = [(s["g0"], s["delta"]) for s in scan_results]
    deltas_scan.sort()
    # Czy delta monotonicznie rosnie/malaje?
    delta_vals = [d for _, d in deltas_scan]
    mono_inc = all(delta_vals[i+1] >= delta_vals[i] - 0.1 for i in range(len(delta_vals)-1))
    mono_dec = all(delta_vals[i+1] <= delta_vals[i] + 0.1 for i in range(len(delta_vals)-1))
    print(f"  delta monotonicznie rosnie: {mono_inc}")
    print(f"  delta monotonicznie maleje: {mono_dec}")
    if not mono_inc and not mono_dec:
        print(f"  delta(g0) jest NIEMONOTONICZNE -- ma strukture.")

# ================================================================
# SECTION 6: Kluczowy test -- Koide K z ekstraktowanych A
# ================================================================
print(f"\n{'='*70}")
print("  6. KOIDE K Z EKSTRAKTOWANYCH A -- weryfikacja")
print("="*70)

if set(results.keys()) >= {"e", "mu", "tau"}:
    # m_i = A_i^4 (cutoff-indep, A_i^4 proporcjonalne do m_i)
    m_e = results["e"]["A4"]
    m_mu = results["mu"]["A4"]
    m_tau = results["tau"]["A4"]

    # Koide K
    s_sum = math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau)
    K = (m_e + m_mu + m_tau) / s_sum**2

    print(f"\n  m_e  (A_e^4)   = {m_e:.8e}")
    print(f"  m_mu (A_mu^4)  = {m_mu:.8e}")
    print(f"  m_tau(A_tau^4) = {m_tau:.8e}")
    print(f"\n  r21 = m_mu/m_e   = {m_mu/m_e:.4f}  (PDG = {M_MU/M_E:.4f})")
    print(f"  r31 = m_tau/m_e  = {m_tau/m_e:.4f}  (PDG = {M_TAU/M_E:.4f})")
    print(f"\n  K = {K:.8f}  (target 2/3 = {2/3:.8f})")
    print(f"  |K - 2/3| = {abs(K - 2/3):.6e}")

    check("T4: Koide K = 2/3 z ekstraktowanych A_i (tol 1%)",
          abs(K - 2/3) < 0.01,
          f"K = {K:.6f}, diff = {abs(K-2/3):.4e}")
    check("T5: r21 zgodne z PDG (tol 1%)",
          abs(m_mu/m_e - M_MU/M_E) / (M_MU/M_E) < 0.01,
          f"r21 = {m_mu/m_e:.4f} vs PDG {M_MU/M_E:.4f}")
    check("T6: r31 zgodne z PDG (tol 5%)",
          abs(m_tau/m_e - M_TAU/M_E) / (M_TAU/M_E) < 0.05,
          f"r31 = {m_tau/m_e:.4f} vs PDG {M_TAU/M_E:.4f}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*70}")
print(f"  TEST REPORT: {PASS} PASS, {FAIL} FAIL out of {PASS+FAIL}")
print(f"{'='*70}")
print(f"""
  INTERPRETACJA:
  - Jesli T1 PASS: delta_i rowniez 120 deg apart -> Z3 na fazach OGONA.
    To oznacza ze Brannen B=sqrt(2) wynikalby z struktury fazowej ogona.
    To byl by silny argument za Koide z ODE.
  - Jesli T1 FAIL: fazy NIE sa 120 deg apart.
    Trzeba szukac innej struktury Z3-rownowaznej.
  - T3 testuje czy A^2(delta) ma forme Brannena.
  - T4 weryfikuje Koide z ekstraktowanych A (kontrola sanity).
""")
