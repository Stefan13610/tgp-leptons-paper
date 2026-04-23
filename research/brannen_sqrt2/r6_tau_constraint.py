#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_tau_constraint.py -- What determines g0^tau?

KEY RESULT from r6_koide_from_ode.py:
  - g0^tau(Koide) = 1.7293, g0^tau/g0^e = 1.989 ~ 2
  - g0_crit = 2.250 (barrier)
  - phi^2*g0^e = 2.276 > g0_crit (BLOCKED)

THREE HYPOTHESES:
  H1: g0^tau = 2*g0^e  (algebraic doubling)
  H2: g0^tau = phi*g0^e + (1-g0^e) = phi*g0^e + delta_e (reflection)
  H3: g0^tau from some ODE-intrinsic property (e.g., energy extremum)

Also explores: is g0^e = 0.86941 itself determined by some condition?
If so, then g0^tau = 2*g0^e would give K = 2/3 as a PREDICTION.

Author: Claudian
Date: 2026-04-15
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
import math

PHI = (1 + math.sqrt(5)) / 2
SQRT2 = math.sqrt(2)

PASS_COUNT = 0
FAIL_COUNT = 0

def check(name, condition, detail=""):
    global PASS_COUNT, FAIL_COUNT
    if condition:
        PASS_COUNT += 1
        print(f"  PASS {name}  {detail}")
    else:
        FAIL_COUNT += 1
        print(f"  FAIL {name}  {detail}")

M_E = 0.510999; M_MU = 105.6584; M_TAU = 1776.86
R21_PDG = M_MU / M_E
R31_PDG = M_TAU / M_E

# ================================================================
# SOLVER (reused)
# ================================================================
def solve_substrate(g0, r_max=300.0, n_points=60000):
    singular = [False]
    def rhs(r, y):
        g, gp = y
        if g < 1e-10:
            singular[0] = True
            g = 1e-10
        if r < 1e-12:
            gpp = (1 - g) / 4.0
        else:
            gpp = (1 - g) - (1.0/g)*gp**2 - (2.0/r)*gp
        return [gp, gpp]
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol, singular[0]

def extract_atail(r, g, r_min=80.0, r_max=250.0):
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 20:
        return None
    u_f = (g[mask] - 1.0) * r_f
    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, u_f, rcond=None)
    A = math.sqrt(coef[0]**2 + coef[1]**2)
    y_hat = coef[0]*np.cos(r_f) + coef[1]*np.sin(r_f)
    rmse = float(np.sqrt(np.mean((u_f - y_hat)**2)))
    return A if rmse/max(A, 1e-10) < 0.05 else None

def get_atail(g0, r_max=300.0):
    sol, sing = solve_substrate(g0, r_max=r_max)
    if sing or not sol.success:
        return None
    return extract_atail(sol.t, sol.y[0])

def get_atail_robust(g0):
    A = get_atail(g0, r_max=300.0)
    if A is not None:
        return A
    for rm, rmin, rmax in [(150, 40, 120), (80, 20, 60), (50, 15, 40)]:
        sol, sing = solve_substrate(g0, r_max=rm, n_points=max(20000, rm*200))
        if not sing and sol.success:
            A = extract_atail(sol.t, sol.y[0], r_min=rmin, r_max=rmax)
            if A is not None:
                return A
    return None

def koide_K(m1, m2, m3):
    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return (m1 + m2 + m3) / s**2

# ================================================================
print("=" * 70)
print("  R6: WHAT DETERMINES g0^tau?")
print("=" * 70)

# Physical values
G0E = 0.86941
G0MU = PHI * G0E
A_E = get_atail(G0E)
A_MU = get_atail(G0MU)
m_e = A_E**4
m_mu = A_MU**4

print(f"\n  g0^e  = {G0E:.5f}")
print(f"  g0^mu = {G0MU:.5f} = phi * g0^e")
print(f"  A_e   = {A_E:.8f}")
print(f"  A_mu  = {A_MU:.8f}")
print(f"  r21   = {(A_MU/A_E)**4:.4f}")

# ================================================================
# SECTION 1: Hypothesis testing - g0^tau candidates
# ================================================================
print(f"\n{'='*70}")
print("  1. g0^tau CANDIDATES")
print("="*70)

g0tau_candidates = {
    "phi^2 * g0e": PHI**2 * G0E,
    "2 * g0e": 2.0 * G0E,
    "phi * g0e + delta_e": PHI * G0E + (1 - G0E),
    "g0e + phi*delta_e": G0E + PHI * (1 - G0E),
    "1 + phi*(g0mu-1)": 1 + PHI * (G0MU - 1),
    "g0mu + delta_e": G0MU + (1 - G0E),
    "g0mu + (g0mu-1)": G0MU + (G0MU - 1),  # = 2*g0mu - 1
    "g0e + 2*delta_e": G0E + 2*(1 - G0E),  # = 2 - g0e
    "2 - g0e": 2.0 - G0E,
    "3 - 2*g0e": 3.0 - 2*G0E,
    "1/g0e + g0e": 1.0/G0E + G0E,
    "phi/g0e": PHI / G0E,
    "phi + delta_e": PHI + (1 - G0E),
}

# The Koide value from previous script
G0TAU_KOIDE = 1.729313

print(f"\n  Koide target: g0^tau = {G0TAU_KOIDE:.6f}")
print(f"\n  {'Candidate':>25s}  {'Value':>10s}  {'Diff':>10s}  {'Diff%':>8s}")
print(f"  {'-'*25}  {'-'*10}  {'-'*10}  {'-'*8}")

for name, val in sorted(g0tau_candidates.items(), key=lambda x: abs(x[1] - G0TAU_KOIDE)):
    diff = val - G0TAU_KOIDE
    diff_pct = abs(diff) / G0TAU_KOIDE * 100
    marker = " <***" if diff_pct < 1 else " <--" if diff_pct < 5 else ""
    print(f"  {name:>25s}  {val:10.6f}  {diff:+10.6f}  {diff_pct:7.3f}%{marker}")

# ================================================================
# SECTION 2: Test each candidate for K = 2/3
# ================================================================
print(f"\n{'='*70}")
print("  2. KOIDE K FOR EACH CANDIDATE")
print("="*70)

print(f"\n  {'Candidate':>25s}  {'g0_tau':>10s}  {'A_tau':>10s}  {'r31':>10s}  {'K':>10s}  {'|K-2/3|':>10s}")
print(f"  {'-'*25}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

best_K_diff = 1.0
best_candidate = ""

for name, val in sorted(g0tau_candidates.items(), key=lambda x: abs(x[1] - G0TAU_KOIDE)):
    if val < 1.1 or val > 3.0:
        continue
    A_tau = get_atail_robust(val)
    if A_tau is not None and A_tau > 0:
        m_tau = A_tau**4
        K = koide_K(m_e, m_mu, m_tau)
        r31 = m_tau / m_e
        K_diff = abs(K - 2/3)
        marker = " <<<" if K_diff < 0.005 else " <--" if K_diff < 0.02 else ""
        print(f"  {name:>25s}  {val:10.5f}  {A_tau:10.6f}  {r31:10.1f}  {K:10.6f}  {K_diff:10.6f}{marker}")
        if K_diff < best_K_diff:
            best_K_diff = K_diff
            best_candidate = name
    else:
        print(f"  {name:>25s}  {val:10.5f}  {'FAILED':>10s}")

print(f"\n  Best candidate: {best_candidate} (|K-2/3| = {best_K_diff:.6f})")


# ================================================================
# SECTION 3: Is g0^e itself determined?
# ================================================================
print(f"\n{'='*70}")
print("  3. WHAT DETERMINES g0^e = 0.86941?")
print("="*70)

print("""
  g0^e = 0.86941 comes from fitting the electron Compton wavelength.
  But is there an ODE-intrinsic property that selects it?

  Possibilities:
  a) g0^e = phi - 1/phi = 0.7236... (NO, doesn't match)
  b) g0^e such that r21 = PDG value (circular)
  c) g0^e from some extremum of the ODE
  d) g0^e from the Compton wavelength formula lambda_C = h/(mc)

  Let's test: does g0^e have any special ODE property?
""")

# Check if g0^e minimizes/maximizes some ODE functional
# Compute energy of soliton as function of g0
def soliton_energy(g0, r_max=100.0, n_pts=20000):
    """Compute E = 4pi * integral[(g^2*g'^2/2 + U(g)-U(1))*r^2 dr]"""
    sol, sing = solve_substrate(g0, r_max=r_max, n_points=n_pts)
    if sing or not sol.success:
        return None
    r = sol.t
    g = sol.y[0]
    gp = sol.y[1]

    # U(g) = g^3/3 - g^4/4;  U(1) = 1/3 - 1/4 = 1/12
    U = g**3/3 - g**4/4
    U1 = 1.0/12
    kinetic = g**2 * gp**2 / 2
    integrand = (kinetic + U - U1) * r**2

    E = 4 * math.pi * np.trapezoid(integrand, r)
    return E

# Also compute A_tail-based mass
print(f"\n  {'g0':>8s}  {'Energy':>12s}  {'A_tail':>10s}  {'A^4':>12s}  {'E/A^4':>10s}")
print(f"  {'-'*8}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*10}")

g0_test = [0.5, 0.6, 0.7, 0.8, 0.85, 0.86941, 0.9, 0.95,
           1.05, 1.1, 1.2, 1.3, 1.4, 1.407, 1.5, 1.6, 1.729, 1.8, 1.9, 2.0]

for g0 in g0_test:
    E = soliton_energy(g0)
    A = get_atail_robust(g0)
    if E is not None and A is not None and A > 1e-10:
        A4 = A**4
        ratio = E / A4 if A4 > 1e-20 else 0
        marker = ""
        if abs(g0 - G0E) < 0.001:
            marker = " <-- electron"
        elif abs(g0 - G0MU) < 0.01:
            marker = " <-- muon"
        elif abs(g0 - G0TAU_KOIDE) < 0.01:
            marker = " <-- tau (Koide)"
        print(f"  {g0:8.4f}  {E:12.6f}  {A:10.6f}  {A4:12.8f}  {ratio:10.4f}{marker}")


# ================================================================
# SECTION 4: The c_M coefficient (E vs A^4)
# ================================================================
print(f"\n{'='*70}")
print("  4. E vs A^4: IS c_M = E/A^4 CONSTANT?")
print("="*70)

print("""
  If m = c_M * A_tail^4, then E/A^4 should give c_M.
  If c_M is CONSTANT across all solitons, mass ~ A^4 is exact.
  If c_M varies, the relationship is approximate.
""")

# Collect E and A^4 for deficit solitons only (E > 0)
deficit_data = []
for g0 in np.linspace(0.3, 0.98, 30):
    E = soliton_energy(g0)
    A = get_atail_robust(g0)
    if E is not None and A is not None and A > 1e-10:
        deficit_data.append((g0, E, A, A**4, E/A**4))

if deficit_data:
    ratios = [d[4] for d in deficit_data]
    r_mean = np.mean(ratios)
    r_std = np.std(ratios)
    r_cv = r_std / r_mean

    print(f"\n  Deficit solitons: c_M = E/A^4")
    print(f"    mean  = {r_mean:.4f}")
    print(f"    std   = {r_std:.4f}")
    print(f"    CV    = {r_cv:.4f} ({r_cv*100:.2f}%)")

    check("T1: c_M approximately constant for deficit solitons",
          r_cv < 0.1,
          f"CV = {r_cv:.4f}")


# ================================================================
# SECTION 5: Deeper exploration of g0^tau = 2*g0^e
# ================================================================
print(f"\n{'='*70}")
print("  5. g0^tau = 2*g0^e HYPOTHESIS")
print("="*70)

print("""
  If g0^tau = 2*g0^e, then:
    delta_e = 1 - g0^e
    delta_mu = phi*g0^e - 1 = (phi-1)*g0^e + g0^e - 1 = phi^(-1)*g0^e + g0^e - 1
    delta_tau = 2*g0^e - 1

  Note: delta_tau = 2*g0^e - 1 = 1 - 2*(1-g0^e) = 1 - 2*delta_e

  The deltas satisfy:
    delta_e + delta_tau = 2*g0^e - 1 + 1 - g0^e = g0^e
    delta_tau = 1 - 2*delta_e

  This gives a SYMMETRIC structure around g0 = 1:
    g0^e  = 1 - delta_e        (deficit)
    g0^tau = 1 + (1 - 2*delta_e)  (excess, at distance 1-2*delta_e from vacuum)
""")

# Test: what g0^e makes K = 2/3 with g0^tau = 2*g0^e?
def K_at_g0e_with_2x_tau(g0e_val):
    """Compute K with g0^mu = phi*g0^e, g0^tau = 2*g0^e"""
    A1 = get_atail(g0e_val)
    A2 = get_atail(PHI * g0e_val)
    A3 = get_atail_robust(2.0 * g0e_val)
    if A1 is None or A2 is None or A3 is None:
        return None
    return koide_K(A1**4, A2**4, A3**4)

print(f"\n  K(g0^e) with g0^tau = 2*g0^e:")
print(f"  {'g0^e':>8s}  {'g0^mu':>8s}  {'g0^tau':>8s}  {'K':>10s}  {'K-2/3':>10s}")

g0e_range = np.linspace(0.50, 0.95, 30)
K_2x = []
g0e_2x = []

for g0e_val in g0e_range:
    K = K_at_g0e_with_2x_tau(g0e_val)
    if K is not None:
        K_2x.append(K)
        g0e_2x.append(g0e_val)
        marker = " <-- physical" if abs(g0e_val - G0E) < 0.01 else ""
        if abs(K - 2/3) < 0.02:
            marker += " <-- K~2/3"
        print(f"  {g0e_val:8.4f}  {PHI*g0e_val:8.4f}  {2*g0e_val:8.4f}  {K:10.6f}  {K-2/3:+10.6f}{marker}")

K_2x = np.array(K_2x)
g0e_2x = np.array(g0e_2x)

# Find crossing of K = 2/3
K_diff_2x = K_2x - 2/3
crossings = np.where(np.diff(np.sign(K_diff_2x)))[0]
if len(crossings) > 0:
    for idx in crossings:
        g0e_cross = g0e_2x[idx] + (2/3 - K_2x[idx]) * (g0e_2x[idx+1] - g0e_2x[idx]) / (K_2x[idx+1] - K_2x[idx])
        print(f"\n  K = 2/3 at g0^e = {g0e_cross:.5f} (with g0^tau = 2*g0^e = {2*g0e_cross:.5f})")
        print(f"  Physical g0^e = {G0E:.5f}")
        print(f"  Difference: {abs(g0e_cross - G0E)/G0E*100:.2f}%")

        check("T2: K=2/3 crossing close to physical g0^e",
              abs(g0e_cross - G0E)/G0E < 0.05,
              f"crossing at {g0e_cross:.5f}, physical at {G0E:.5f}")
else:
    print(f"\n  No K = 2/3 crossing found in range!")


# ================================================================
# SECTION 6: Alternative ladder: e*phi->mu, mu*(2/phi)->tau
# ================================================================
print(f"\n{'='*70}")
print("  6. ALTERNATIVE LADDERS")
print("="*70)

print("""
  Instead of phi-ladder for ALL steps, what if:
    g0^mu = phi * g0^e     (confirmed by r21)
    g0^tau = R * g0^mu     (unknown R for tau step)

  R = g0^tau / g0^mu = 1.729 / 1.407 = 1.229
  What is R algebraically?
""")

R_tau_mu = G0TAU_KOIDE / G0MU
print(f"  R = g0^tau(Koide)/g0^mu = {R_tau_mu:.6f}")

# Check algebraic values
candidates_R = [
    (1.0, "1"),
    (SQRT2 - 0.2, "sqrt(2)-0.2"),
    (5.0/4, "5/4"),
    (PHI - 1.0/3, "phi-1/3"),
    (1 + 1.0/PHI**2, "1+1/phi^2"),
    (PHI/math.sqrt(PHI+1), "phi/sqrt(phi+1)"),
    (math.sqrt(PHI), "sqrt(phi)"),
    (PHI**2 / 2, "phi^2/2"),
    (2.0/PHI, "2/phi"),
    (3.0/PHI - 1, "3/phi-1"),
    (math.sqrt(3)/math.sqrt(2), "sqrt(3/2)"),
    (1 + 1.0/(PHI+1), "1+1/(phi+1)"),
    (PHI/(1+1.0/PHI), "phi/(1+1/phi)"),
    ((PHI+1)/math.sqrt(PHI**2+PHI+1), "(phi+1)/sqrt(phi^2+phi+1)"),
    (math.log(PHI) + 1, "1+ln(phi)"),
]

print(f"\n  Closest algebraic values to R = {R_tau_mu:.6f}:")
for val, name in sorted(candidates_R, key=lambda x: abs(x[0] - R_tau_mu)):
    diff = abs(val - R_tau_mu)
    if diff < 0.1:
        print(f"    {name:>25s} = {val:.6f}  (diff = {diff:.6f})")


# ================================================================
# SECTION 7: The g0^e self-consistency condition
# ================================================================
print(f"\n{'='*70}")
print("  7. SELF-CONSISTENCY: g0^e FROM COMPTON WAVELENGTH")
print("="*70)

print("""
  In TGP, g0^e is determined by the electron Compton wavelength:
    lambda_C = h / (m_e * c)

  The soliton radius R_sol is set by the ODE length scale.
  The physical mass is m = c_M * A_tail^4.

  So lambda_C(g0^e) = h / (c_M * A_tail(g0^e)^4 * c)

  The constraint is: lambda_C = R_sol (soliton = Compton wavelength).
  This determines g0^e.

  But this requires knowing c_M and the length scale.
  In DIMENSIONLESS terms: the ratio r21 = (A_mu/A_e)^4 is independent
  of c_M and the length scale. So r21 is a PURE number predicted by
  the ODE + phi-ladder.

  KEY QUESTION: Does the ODE predict r21 = 206.768 for ANY g0^e,
  or only for g0^e = 0.86941?
""")

# Compute r21 as function of g0^e
print(f"\n  r21(g0^e) = [A_tail(phi*g0^e) / A_tail(g0^e)]^4:")
print(f"  {'g0^e':>8s}  {'A_e':>10s}  {'A_mu':>10s}  {'r21':>10s}  {'r21/PDG':>8s}")

r21_data = []
for g0e_val in np.linspace(0.30, 0.95, 40):
    g0mu_val = PHI * g0e_val
    A1 = get_atail(g0e_val)
    A2 = get_atail(g0mu_val)
    if A1 is not None and A2 is not None and A1 > 1e-10:
        r21_val = (A2/A1)**4
        r21_data.append((g0e_val, A1, A2, r21_val))
        marker = " <-- physical" if abs(g0e_val - G0E) < 0.01 else ""
        if g0e_val in [0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95] or abs(g0e_val - G0E) < 0.02:
            print(f"  {g0e_val:8.4f}  {A1:10.6f}  {A2:10.6f}  {r21_val:10.2f}  {r21_val/R21_PDG:8.4f}{marker}")

if r21_data:
    r21_vals = [d[3] for d in r21_data]
    r21_mean = np.mean(r21_vals)
    r21_std = np.std(r21_vals)
    r21_cv = r21_std / r21_mean

    print(f"\n  r21 statistics:")
    print(f"    mean = {r21_mean:.2f}")
    print(f"    std  = {r21_std:.2f}")
    print(f"    CV   = {r21_cv:.4f}")

    if r21_cv > 0.05:
        print(f"    r21 VARIES strongly with g0^e (NOT universal!)")
        print(f"    r21 at physical g0^e is a consequence of g0^e = 0.86941")

        # Find g0^e that gives PDG r21
        g0e_arr = np.array([d[0] for d in r21_data])
        r21_arr = np.array([d[3] for d in r21_data])
        r21_diff = r21_arr - R21_PDG
        cross = np.where(np.diff(np.sign(r21_diff)))[0]
        if cross.size > 0:
            idx = cross[0]
            g0e_pdg = g0e_arr[idx] + (R21_PDG - r21_arr[idx]) * (g0e_arr[idx+1] - g0e_arr[idx]) / (r21_arr[idx+1] - r21_arr[idx])
            print(f"\n    r21 = PDG at g0^e = {g0e_pdg:.5f}")
            print(f"    Actual g0^e = {G0E:.5f}")
            print(f"    Agreement: {abs(g0e_pdg - G0E)/G0E*100:.2f}%")
    else:
        print(f"    r21 approximately CONSTANT (universal!)")


# ================================================================
# SECTION 8: Summary chain
# ================================================================
print(f"\n{'='*70}")
print("  8. THE FULL DERIVATION CHAIN (STATUS)")
print("="*70)

print(f"""
  WHAT WE CAN DERIVE:
  1. g_ij = g*delta_ij  ->  substrate ODE (alpha=1)       [PROVEN]
  2. g0^e = 0.86941     <-  Compton wavelength matching    [INPUT: lambda_C]
  3. g0^mu = phi*g0^e   <-  phi-ladder                     [ASSUMPTION: phi spacing]
  4. r21 = 206.55       <-  ODE + steps 2+3                [DERIVED: 0.10% from PDG]
  5. g0^tau = 1.729     <-  ???                             [OPEN]
  6. K = 2/3            <-  depends on g0^tau              [OPEN: equivalent to step 5]
  7. B = sqrt(2)        <-  algebraic from K=2/3           [PROVEN if K=2/3]
  8. N = 3              <-  barrier + Koide g0^tau < g0_crit [PROVEN if K=2/3]

  THE GAP: Step 5. What determines g0^tau?

  CANDIDATES:
  A. g0^tau = 2*g0^e              (factor-2 ladder)  ->  K differs from 2/3 by ~?
  B. g0^tau = g0_crit             (barrier)          ->  K = ?
  C. g0^tau from extremum of some functional          ->  needs investigation
  D. g0^tau = phi*g0^e + delta_e  (shifted ladder)    ->  K = ?
  E. K = 2/3 as INPUT (from Z3 symmetry, GL(3,F2))   ->  g0^tau = 1.729

  If option E is correct, then K=2/3 is a SYMMETRY PRINCIPLE,
  and the full chain requires TWO inputs:
    INPUT 1: g0^e (from lambda_C = soliton radius)
    INPUT 2: K = 2/3 (from Z3 symmetry)

  The phi-ladder gives ONE mass ratio (r21).
  Koide gives the SECOND (r31).
  Together: all 3 masses from g0^e + phi + K.
""")

print(f"\n{'='*70}")
print(f"  TEST REPORT: {PASS_COUNT} PASS, {FAIL_COUNT} FAIL out of {PASS_COUNT + FAIL_COUNT}")
print(f"{'='*70}")
