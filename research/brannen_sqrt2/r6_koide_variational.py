#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_koide_variational.py -- Czy K=2/3 jest wartoscia stacjonarna jakiegos funkcjonalu?

Po negatywnym wyniku r6_tail_phase_z3.py (fazy NIE sa Z3-symetryczne),
nowa hipoteza: istnieje funkcjonal S[g0^e, g0^mu, g0^tau] ktory:
  1. Osiaga ekstremum lub jest stacjonarny przy fizycznych wartosciach
  2. Wymusza K=2/3 jako konsekwencje

POTENCJALNE FUNKCJONALY:
  S1 = suma energii solitonowych E(g0^e) + E(g0^mu) + E(g0^tau)
  S2 = log(det) = ln(m_e * m_mu * m_tau) -- wyznacznik masowy
  S3 = entropia von Neumanna rozkladu mas: -Σ p_i ln p_i  z p_i = m_i/Σm_i
  S4 = sum(A_i * sqrt(m_i)) -- Peters-Young type
  S5 = Koide Q_K wprost -- trywialne
  S6 = funkcjonal typu "free energy" z constraint na masy leptonow

HIPOTEZA: Fizyczne g0^e, g0^mu, g0^tau SA punktami stacjonarnymi S
z wiezami:
  - g0^mu = phi * g0^e (phi-drabinka)
  - g0^e = ustalone przez Compton

To redukuje problem do 1D: czy S(g0^tau | g0^e, g0^mu) ma ekstremum przy 1.729?

Drugi atak: CV(sqrt(m)) = 1 (rownowazne K=2/3).
Czy CV=1 jest ekstremum funkcjonalu, na przyklad dla sztywnie wybranych g0^e, g0^mu,
CV(g0^tau) ma minimum/maximum przy 1.729?

Author: Claudian
Date: 2026-04-16
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq
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
    print(f"  {mark}  {name:45s}  {detail}")

G0_E   = 0.86941
G0_MU  = PHI * G0_E
G0_TAU_KOIDE = 1.7293  # z inwersji Koide K=2/3

M_E, M_MU, M_TAU = 0.510999, 105.6584, 1776.86

def rhs(r, y):
    g, gp = y
    if g < 1e-10: g = 1e-10
    if r < 1e-12:
        gpp = (1 - g) / 4.0
    else:
        gpp = (1 - g) - (1.0/g)*gp**2 - (2.0/r)*gp
    return [gp, gpp]

def solve_substrate(g0, r_max=300.0, n_points=60000):
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol

def fit_A(r, g, r_min=80.0, r_max=250.0):
    mask = (r >= r_min) & (r <= r_max)
    if mask.sum() < 20: return None
    r_f, u_f = r[mask], (g[mask]-1.0)*r[mask]
    X = np.column_stack([np.sin(r_f), np.cos(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, u_f, rcond=None)
    A = math.sqrt(coef[0]**2 + coef[1]**2)
    y_hat = coef[0]*np.sin(r_f) + coef[1]*np.cos(r_f)
    if np.sqrt(np.mean((u_f-y_hat)**2))/max(A,1e-12) > 0.05: return None
    return A

def get_A(g0):
    sol = solve_substrate(g0)
    if not sol.success: return None
    return fit_A(sol.t, sol.y[0])

def soliton_energy(g0):
    sol = solve_substrate(g0, r_max=100.0, n_points=20000)
    if not sol.success: return None
    r = sol.t; g = sol.y[0]; gp = sol.y[1]
    U = g**3/3 - g**4/4
    U1 = 1.0/12
    kinetic = g**2 * gp**2 / 2
    integrand = (kinetic + U - U1) * r**2
    return 4 * math.pi * np.trapezoid(integrand, r)

def koide_K(m_e, m_mu, m_tau):
    s = math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau)
    return (m_e + m_mu + m_tau) / s**2

def CV_sqrt_m(m_e, m_mu, m_tau):
    sq = np.array([math.sqrt(m_e), math.sqrt(m_mu), math.sqrt(m_tau)])
    mu = sq.mean()
    sd = sq.std()  # populacja sigma
    return sd / mu

# ================================================================
print("=" * 70)
print("  R6: CZY K=2/3 JEST STACJONARNA WARTOSCIA FUNKCJONALU?")
print("=" * 70)

A_E = get_A(G0_E)
A_MU = get_A(G0_MU)
print(f"\n  g0^e  = {G0_E:.5f}, A_e  = {A_E:.6f}")
print(f"  g0^mu = {G0_MU:.5f}, A_mu = {A_MU:.6f}")
print(f"  g0^tau_Koide = {G0_TAU_KOIDE:.5f}")

# ================================================================
# SECTION 1: K(g0^tau) jako funkcja g0^tau
# ================================================================
print(f"\n{'='*70}")
print("  1. SKAN K(g0^tau) PRZY USTALONYCH g0^e, g0^mu")
print("="*70)

g0tau_scan = np.linspace(1.40, 2.15, 80)
K_vals = []
A_tau_vals = []
E_tau_vals = []
CV_vals = []
g0_keep = []

for g0 in g0tau_scan:
    A = get_A(g0)
    if A is None or A <= 0:
        continue
    m_e, m_mu, m_tau = A_E**4, A_MU**4, A**4
    K = koide_K(m_e, m_mu, m_tau)
    CV = CV_sqrt_m(m_e, m_mu, m_tau)
    E = soliton_energy(g0)
    K_vals.append(K)
    A_tau_vals.append(A)
    E_tau_vals.append(E)
    CV_vals.append(CV)
    g0_keep.append(g0)

g0_keep = np.array(g0_keep)
K_vals = np.array(K_vals)
CV_vals = np.array(CV_vals)
A_tau_vals = np.array(A_tau_vals)
E_tau_vals = np.array([e if e is not None else np.nan for e in E_tau_vals])

print(f"\n  Sample {'g0^tau':>8s} {'A_tau':>10s} {'m_tau/m_e':>11s} {'K':>10s} {'K-2/3':>10s} {'CV(sqm)':>10s} {'E_tau':>12s}")
print(f"  ------ {'-'*8} {'-'*10} {'-'*11} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
for i in range(0, len(g0_keep), 5):
    g0 = g0_keep[i]
    marker = ""
    if abs(g0 - G0_TAU_KOIDE) < 0.02: marker = " <-- Koide"
    print(f"  [{i:3d}]  {g0:8.4f} {A_tau_vals[i]:10.6f} {A_tau_vals[i]**4/A_E**4:11.2f} {K_vals[i]:10.6f} {K_vals[i]-2/3:+10.6f} {CV_vals[i]:10.6f} {E_tau_vals[i]:12.6f}{marker}")

# Znajdz minimum |K - 2/3|
idx_best = int(np.argmin(np.abs(K_vals - 2/3)))
g0_K23 = g0_keep[idx_best]
K_best = K_vals[idx_best]
print(f"\n  Minimum |K - 2/3|: g0^tau = {g0_K23:.5f}, K = {K_best:.8f}")

# Znajdz minimum |CV - 1| -- tozsame z K=2/3
idx_cv = int(np.argmin(np.abs(CV_vals - 1.0)))
g0_CV1 = g0_keep[idx_cv]
CV_best = CV_vals[idx_cv]
print(f"  Minimum |CV - 1|   : g0^tau = {g0_CV1:.5f}, CV = {CV_best:.8f}")
print(f"\n  Target (z inwersji Koide): g0^tau = {G0_TAU_KOIDE:.5f}")

check("T1: g0_K23 pasuje do Koide inversion (tol 1%)",
      abs(g0_K23 - G0_TAU_KOIDE)/G0_TAU_KOIDE < 0.01,
      f"g0_K23 = {g0_K23:.5f}")

# ================================================================
# SECTION 2: ekstrema innych funkcjonalow
# ================================================================
print(f"\n{'='*70}")
print("  2. EKSTREMUM INNYCH FUNKCJONALOW")
print("="*70)

# log-determinant masowy
log_det = np.log(A_E**4 * A_MU**4 * A_tau_vals**4)
# Peters-Young type: sum(A_i * sqrt(m_i)) = sum(A_i * A_i^2) = sum(A^3)
PY = A_E**3 + A_MU**3 + A_tau_vals**3
# Entropia Shannona mas
ms = np.array([(A_E**4, A_MU**4, at**4) for at in A_tau_vals])
p = ms / ms.sum(axis=1, keepdims=True)
H_shannon = -np.sum(p * np.log(p + 1e-20), axis=1)

# Suma energii
E_sum = soliton_energy(G0_E) + soliton_energy(G0_MU) + E_tau_vals

# Wariancja log mass
log_ms = np.log(ms)
var_log = log_ms.var(axis=1)

print(f"\n  Ekstrema funkcjonalow wzgledem g0^tau:")
print(f"  {'Funkcjonal':>20s}  {'g0^tau ekstrem':>14s}  {'Wartosc':>12s}  {'g0_Koide diff':>14s}")
print(f"  {'-'*20}  {'-'*14}  {'-'*12}  {'-'*14}")

candidates = {
    "log(det m)": (log_det, False),  # mono rosnie, brak ekstremum
    "sum(A^3)":   (PY, False),
    "H Shannon":  (H_shannon, True),   # ma maksimum
    "E_sum":      (E_sum, True),
    "var(log m)": (var_log, True),     # ma minimum/maksimum
    "K":          (K_vals, True),      # K ma maksimum!
    "CV":         (CV_vals, True),
}

for name, (vals, find_ext) in candidates.items():
    vals = np.array(vals)
    if not find_ext:
        continue
    # Znajdz lokalny max i min
    i_max = int(np.argmax(vals))
    i_min = int(np.argmin(vals))
    g_max = g0_keep[i_max]; v_max = vals[i_max]
    g_min = g0_keep[i_min]; v_min = vals[i_min]

    # Ktory jest blizszy Koide?
    diff_max = abs(g_max - G0_TAU_KOIDE)
    diff_min = abs(g_min - G0_TAU_KOIDE)
    if diff_max < diff_min:
        label = "MAX"
        g_ext = g_max
        v_ext = v_max
        diff = diff_max
    else:
        label = "MIN"
        g_ext = g_min
        v_ext = v_min
        diff = diff_min
    marker = " <--" if diff < 0.05 else ""
    print(f"  {name:>20s}  {g_ext:10.5f} ({label}) {v_ext:12.6f}  {diff:10.6f}{marker}")

# K ma maximum = 2/3 przy g0_K23 -- to jest tautologia poniewaz K=2/3 jest max mozliwy przy Koide
# Ale geometrycznie K ogranicza sie od gory przez 1 i od dolu przez 1/3

# ================================================================
# SECTION 3: K CONSTRAINED by CAUCHY-SCHWARTZ - czy fizyczny stan to MAX K?
# ================================================================
print(f"\n{'='*70}")
print("  3. KAUCHY-SCHWARZ: CZY K MIESCI SIE W [1/3, 1]?")
print("="*70)

print("""
  Z nierownosci Cauchy'ego-Schwartza:
    (sum sqrt(m_i))^2 <= N * sum(m_i)    =>   K >= 1/N = 1/3
    K = 1 iff wszystkie m_i rowne (zdegenerowane)
    K = 1/N iff jedno m_i dominuje (1 generacja)

  Maksimum K to 1, minimum to 1/N = 1/3.
  Dlaczego K = 2/N = 2/3 dokladnie? 2/3 = SREDNIA arytmetyczna (1/3, 1).

  HIPOTEZA: K = 2/3 jest PUNKT POSREDNI miedzy maksymalna hierarchia (1/3)
  a pelna degeneracja (1). Jakas zasada "najmniej narzuconego" wybor?
""")

# Skan K dla wszystkich 3 leptonow jako funkcja 1 g0 przy ustalonych 2 pozostalych
# Pokaz zakres K
print(f"\n  Zakres K w skanie g0^tau:")
print(f"    min K = {K_vals.min():.6f} przy g0^tau = {g0_keep[K_vals.argmin()]:.5f}")
print(f"    max K = {K_vals.max():.6f} przy g0^tau = {g0_keep[K_vals.argmax()]:.5f}")
print(f"    K Koide target = 2/3 = {2/3:.6f}")
print(f"    K przy g0^tau_Koide = {K_vals[np.argmin(np.abs(g0_keep - G0_TAU_KOIDE))]:.6f}")

# Dwa punkty gdzie K przekracza 2/3 -- to jest CROSS
crossings = []
for i in range(len(K_vals)-1):
    if (K_vals[i] - 2/3)*(K_vals[i+1] - 2/3) < 0:
        # Interpolacja
        frac = (2/3 - K_vals[i]) / (K_vals[i+1] - K_vals[i])
        g_cross = g0_keep[i] + frac * (g0_keep[i+1] - g0_keep[i])
        crossings.append(g_cross)

print(f"\n  Liczba przeciec K = 2/3: {len(crossings)}")
for j, g in enumerate(crossings):
    print(f"    Przeciecie {j+1}: g0^tau = {g:.5f}  (diff od Koide = {g - G0_TAU_KOIDE:+.5f})")

# ================================================================
# SECTION 4: K jako funkcja g0_e (zmienna E) z g0_mu = phi*g0_e, g0_tau = 2*g0_e
# ================================================================
print(f"\n{'='*70}")
print("  4. K JAKO F(g0^e) Z g0^mu = phi*g0^e, g0^tau = Koide_inversion")
print("="*70)

print("""
  Hipoteza: Dla KAZDEGO g0^e (w zakresie dopuszczalnym) istnieje g0^tau taki ze K=2/3.
  Kluczowe pytanie: czy ten g0^tau ma struktura algebraiczna vs g0^e?
""")

def find_g0_tau_for_K23(g0_e):
    """Znajdz g0^tau takie ze K(g0^e, phi*g0^e, g0^tau) = 2/3"""
    A_e = get_A(g0_e)
    A_mu = get_A(PHI * g0_e)
    if A_e is None or A_mu is None:
        return None, None, None

    def fun(g0_tau):
        A_t = get_A(g0_tau)
        if A_t is None or A_t <= 0:
            return 1.0
        return koide_K(A_e**4, A_mu**4, A_t**4) - 2/3

    # Znajdz roots w roznych zakresach g0_tau
    try:
        g0_tau = brentq(fun, 1.2, 2.1, xtol=1e-5)
        A_t = get_A(g0_tau)
        return g0_tau, A_t, A_t**4 / A_e**4 if A_e**4 > 0 else None
    except Exception:
        return None, None, None

print(f"\n  {'g0^e':>8s} {'g0^mu':>8s} {'g0^tau':>8s} {'r31':>10s} {'g0^tau/g0^e':>12s} {'g0^tau/g0^mu':>14s}")
print(f"  {'-'*8} {'-'*8} {'-'*8} {'-'*10} {'-'*12} {'-'*14}")

g0_e_scan = [0.5, 0.6, 0.7, 0.8, G0_E, 0.9, 0.95]
ratios = []
for g0e in g0_e_scan:
    g0_mu = PHI * g0e
    if g0_mu > 2.15:
        continue
    g0t, A_t, r31 = find_g0_tau_for_K23(g0e)
    if g0t is None:
        print(f"  {g0e:8.4f} {g0_mu:8.4f}  FAILED")
        continue
    r_e = g0t / g0e
    r_mu = g0t / g0_mu
    ratios.append((g0e, g0t, r_e, r_mu))
    marker = " <-- physical" if abs(g0e - G0_E) < 0.01 else ""
    print(f"  {g0e:8.4f} {g0_mu:8.4f} {g0t:8.4f} {r31:10.2f} {r_e:12.6f} {r_mu:14.6f}{marker}")

# Stale stosunki?
if ratios:
    r_e_vals = np.array([r[2] for r in ratios])
    r_mu_vals = np.array([r[3] for r in ratios])
    print(f"\n  Stosunki g0^tau/g0^e:")
    print(f"    mean = {r_e_vals.mean():.6f}, std = {r_e_vals.std():.6f}, CV = {r_e_vals.std()/r_e_vals.mean():.4f}")
    print(f"  Stosunki g0^tau/g0^mu:")
    print(f"    mean = {r_mu_vals.mean():.6f}, std = {r_mu_vals.std():.6f}, CV = {r_mu_vals.std()/r_mu_vals.mean():.4f}")

    check("T2: g0^tau/g0^e prawie stale (CV<5%)",
          r_e_vals.std()/r_e_vals.mean() < 0.05,
          f"CV = {r_e_vals.std()/r_e_vals.mean():.4f}")
    check("T3: g0^tau/g0^mu prawie stale (CV<5%)",
          r_mu_vals.std()/r_mu_vals.mean() < 0.05,
          f"CV = {r_mu_vals.std()/r_mu_vals.mean():.4f}")

# ================================================================
# SECTION 5: A relacje kwadratowe -- czy A_mu^2 = A_e * A_tau lub podobne?
# ================================================================
print(f"\n{'='*70}")
print("  5. RELACJE KWADRATOWE W A_i")
print("="*70)

A_tau_koide = get_A(G0_TAU_KOIDE)
if A_tau_koide is not None:
    print(f"\n  A_e   = {A_E:.8f}")
    print(f"  A_mu  = {A_MU:.8f}")
    print(f"  A_tau = {A_tau_koide:.8f}")

    # Test: A_mu^2 = A_e * A_tau? (geometric progression)
    amu_sq = A_MU**2
    ae_times_atau = A_E * A_tau_koide
    print(f"\n  A_mu^2       = {amu_sq:.8f}")
    print(f"  A_e * A_tau  = {ae_times_atau:.8f}")
    print(f"  ratio        = {amu_sq/ae_times_atau:.6f}")

    # A_tau - A_mu vs A_mu - A_e (arithmetic progression?)
    d1 = A_MU - A_E
    d2 = A_tau_koide - A_MU
    print(f"\n  A_mu - A_e   = {d1:.8f}")
    print(f"  A_tau - A_mu = {d2:.8f}")
    print(f"  ratio        = {d2/d1:.6f}")

    # A_e + A_tau vs 2*A_mu (arithmetic mean?)
    print(f"\n  A_e + A_tau  = {A_E + A_tau_koide:.8f}")
    print(f"  2 * A_mu     = {2*A_MU:.8f}")
    print(f"  ratio        = {(A_E + A_tau_koide)/(2*A_MU):.6f}")

    # sqrt(m_e) + sqrt(m_tau) vs 2*sqrt(m_mu) -- to jest Brannen!
    sme = A_E**2; smmu = A_MU**2; smtau = A_tau_koide**2
    print(f"\n  sqrt(m_e) = A_e^2   = {sme:.8f}")
    print(f"  sqrt(m_mu) = A_mu^2  = {smmu:.8f}")
    print(f"  sqrt(m_tau) = A_tau^2= {smtau:.8f}")

    # Brannen: sqrt(m_i) = a(1 + sqrt(2)*cos(theta + 2pi*i/3))
    # Implikuje: sum sqrt(m) = 3a (bo sum cos = 0)
    # Implikuje: sqrt(m_e) + sqrt(m_tau) + sqrt(m_mu) = 3a
    # oraz: sqrt(m_e)*sqrt(m_tau) + sqrt(m_mu)*sqrt(m_tau) + sqrt(m_e)*sqrt(m_mu) = 3a^2 - (3a^2)/2 = 3a^2/2 (od TOZSAMOSC cos)
    # Pokazmy: czy a = (sum sqrt(m))/3 dobrze przewiduje Brannen?
    a_brannen = (sme + smmu + smtau) / 3

    # Czysto: s_e*s_mu + s_e*s_tau + s_mu*s_tau
    s_e, s_mu, s_tau = sme, smmu, smtau
    sum_pairs = s_e*s_mu + s_e*s_tau + s_mu*s_tau
    # Brannen TOZSAMOSC: sum_pairs = 3a^2 / 2  (bo cos^2 sum = 3/2)
    # Rigorously: sum (cos_i + cos_j)^2 = ...; see r6_fourier_z3_proof
    #
    # Najczystsze: B^2 = 2 <=> sum sqrt(m)^2 = 2 * (sum sqrt(m))^2 / 3
    # tozsame z K = 2/3: Sum sqrt(m)^2 = Sum m; so to jest po prostu K = Sum m / (Sum sqrt m)^2 = 2/3
    K_check = (s_e**2 + s_mu**2 + s_tau**2) / (s_e + s_mu + s_tau)**2
    print(f"\n  Sprawdzenie Koide K z A:")
    print(f"    K = (A_e^4+A_mu^4+A_tau^4)/(A_e^2+A_mu^2+A_tau^2)^2 = {K_check:.8f}")
    print(f"    target 2/3 = {2/3:.8f}")
    print(f"    diff = {abs(K_check - 2/3):.4e}")

    check("T4: Koide K=2/3 z A przy g0_Koide (tol 1%)",
          abs(K_check - 2/3) < 0.01,
          f"K = {K_check:.6f}")

# ================================================================
# SECTION 6: Klucz: czy istnieje funkcjonal S majacy EKSTREMUM dokladnie przy K=2/3?
# ================================================================
print(f"\n{'='*70}")
print("  6. POSZUKIWANIE EKSTREMUM S[g0^tau]")
print("="*70)

print("""
  Buduje rozne kandydaty S[g0^tau] i sprawdzam czy maja ekstremum przy g0^tau=1.729.

  Krytycznie: S musi miec UZASADNIENIE FIZYCZNE (nie wyglad retrofit).

  Kandydaci fizyczni:
  a) Energia solitonowa suma: E_e + E_mu + E_tau
  b) L2 norma masy: sqrt(Sum m_i^2) -- "total mass energy"
  c) Free energy type: -log(Z) gdzie Z = sum exp(-m_i)
  d) Fisher info: sum (sqrt(m_i+1) - sqrt(m_i))^2
  e) Wigner-Yanase-Dyson skew-information dla rozkladu mas
""")

# Sprawdz czy dowolny z tych ma ekstremum przy g0_tau_Koide
print(f"\n  {'Funkcjonal':>25s}  {'g0 ekstremum':>14s}  {'Wartosc':>12s}  {'Typ':>5s}  {'diff Koide':>11s}")
print(f"  {'-'*25}  {'-'*14}  {'-'*12}  {'-'*5}  {'-'*11}")

m_arr = A_tau_vals**4
me_fix = A_E**4
mmu_fix = A_MU**4

funcs = {
    "E_e+E_mu+E_tau":     soliton_energy(G0_E) + soliton_energy(G0_MU) + E_tau_vals,
    "sum m_i^2 (L2)":     me_fix**2 + mmu_fix**2 + m_arr**2,
    "-log(sum exp(-m))":  -np.log(np.exp(-me_fix) + np.exp(-mmu_fix) + np.exp(-m_arr)),
    "Fisher info":        (np.sqrt(me_fix+1)-np.sqrt(me_fix))**2 + (np.sqrt(mmu_fix+1)-np.sqrt(mmu_fix))**2 + (np.sqrt(m_arr+1)-np.sqrt(m_arr))**2,
    "Sum log m":          np.log(me_fix) + np.log(mmu_fix) + np.log(m_arr),
    "sigma_log(m)^2":     np.array([np.var(np.log([me_fix, mmu_fix, m])) for m in m_arr]),
    "sum sqrt(m)":        np.sqrt(me_fix) + np.sqrt(mmu_fix) + np.sqrt(m_arr),
    "sum m^(2/3)":        me_fix**(2/3) + mmu_fix**(2/3) + m_arr**(2/3),
    "m_tau/m_mu * m_mu/m_e": (m_arr/mmu_fix) * (mmu_fix/me_fix),  # r31
    "K_Koide":            K_vals,
}

for name, vals in funcs.items():
    vals = np.array(vals, dtype=float)
    i_max = int(np.argmax(vals))
    i_min = int(np.argmin(vals))
    g_max = g0_keep[i_max]; g_min = g0_keep[i_min]
    diff_max = abs(g_max - G0_TAU_KOIDE)
    diff_min = abs(g_min - G0_TAU_KOIDE)
    # Szukamy EKSTREMUM (nie brzegowego)
    # Sprawdz: czy g_max i g_min sa punktami wewnetrznymi?
    is_max_interior = 0 < i_max < len(vals) - 1
    is_min_interior = 0 < i_min < len(vals) - 1

    if is_max_interior and diff_max < 0.05:
        print(f"  {name:>25s}  {g_max:14.5f}  {vals[i_max]:12.6f}  {'MAX':>5s}  {diff_max:11.6f}  <-- !!")
    elif is_min_interior and diff_min < 0.05:
        print(f"  {name:>25s}  {g_min:14.5f}  {vals[i_min]:12.6f}  {'MIN':>5s}  {diff_min:11.6f}  <-- !!")
    elif is_max_interior:
        print(f"  {name:>25s}  {g_max:14.5f}  {vals[i_max]:12.6f}  {'MAX':>5s}  {diff_max:11.6f}")
    elif is_min_interior:
        print(f"  {name:>25s}  {g_min:14.5f}  {vals[i_min]:12.6f}  {'MIN':>5s}  {diff_min:11.6f}")
    else:
        print(f"  {name:>25s}  (brzeg)         {'':12s}  {'':>5s}  {'':>11s}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*70}")
print(f"  TEST REPORT: {PASS} PASS, {FAIL} FAIL out of {PASS+FAIL}")
print(f"{'='*70}")
