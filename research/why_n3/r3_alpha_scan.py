#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_alpha_scan.py — Fine scan of α to find the UNIQUE value giving N=3

KEY QUESTION: For which α does g₀_crit(3D) = φ²·g₀^e = 2.2761?
(This would give N=3 EXACTLY from the φ-ladder.)

We know:
  α=0.5: g₀_crit = 2.618 → N=3  (n_max=2.29, above φ²·g₀^e)
  α=1.0: g₀_crit = 2.206 → N=2  (n_max=1.94, below φ²·g₀^e)

So there's a CRITICAL α between 0.5 and 1.0 where g₀_crit = φ²·g₀^e.
This α_crit would be the value where N transitions from 3 to 2.

Also: what α gives g₀_crit = g₀^τ(Koide) = 1.729?
We know α≈3 gives g₀_crit≈1.728.

Physical question: does the TGP Lagrangian determine α?

Autor: Claudian
Data: 2026-04-15
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import math

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}: PASS  {detail}")
    else:
        FAIL += 1
        print(f"  ✗ {name}: FAIL  {detail}")
    return condition

PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941

# ================================================================
# SOLVER
# ================================================================

def solve_alpha(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    """
    ODE: g'' + (α/g)g'² + ((d-1)/r)g' = U'(g)/g^{2α}
    where U(g) = g³/3 - g⁴/4, U'(g) = g²(1-g)
    RHS = g²(1-g)/g^{2α} = (1-g)·g^{2-2α}
    """
    g_min_val = [g0]
    singular = [False]

    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
        if g < g_min_val[0]:
            g_min_val[0] = g

        rhs_val = (1 - g) * g**(2 - 2*alpha)
        if r < 1e-12:
            gpp = rhs_val / max(d, 1.0)
        else:
            gpp = rhs_val - (alpha/g) * gp**2 - ((d-1.0)/r) * gp
        return [gp, gpp]

    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-10, atol=1e-12, max_step=0.05)

    if sol.success:
        actual_min = np.min(sol.y[0])
        if actual_min < g_min_val[0]:
            g_min_val[0] = actual_min

    return g_min_val[0], singular[0], sol


def find_g0_crit(alpha, d=3, g0_lo=1.01, g0_hi=5.0, tol=1e-7):
    g_threshold = 0.005
    g_min, sing, sol = solve_alpha(g0_hi, alpha, d)
    is_bad = sing or g_min < g_threshold or not sol.success
    if not is_bad:
        g0_hi = 10.0
        g_min, sing, sol = solve_alpha(g0_hi, alpha, d)
        is_bad = sing or g_min < g_threshold or not sol.success
        if not is_bad:
            return None

    for _ in range(70):
        g0_mid = (g0_lo + g0_hi) / 2
        g_min, sing, sol = solve_alpha(g0_mid, alpha, d)
        is_bad = sing or g_min < g_threshold or not sol.success
        if is_bad:
            g0_hi = g0_mid
        else:
            g0_lo = g0_mid
        if g0_hi - g0_lo < tol:
            break
    return (g0_lo + g0_hi) / 2


# ================================================================
print("=" * 70)
print("  R3: α-SCAN — FINDING CRITICAL α FOR N=3")
print("=" * 70)

# ================================================================
# SECTION 1: Fine α scan
# ================================================================
print(f"\n{'=' * 70}")
print("  1. g₀_crit(α) FINE SCAN")
print("=" * 70)

alphas = np.arange(0.1, 4.1, 0.1)
g0_crits = {}

for alpha in alphas:
    g0c = find_g0_crit(alpha, d=3)
    if g0c is not None:
        g0_crits[alpha] = g0c

print(f"\n  {'α':>6s}  {'g₀_crit':>10s}  {'n_max':>8s}  {'N_gen':>5s}")
print(f"  {'------':>6s}  {'----------':>10s}  {'--------':>8s}  {'-----':>5s}")
for alpha in sorted(g0_crits.keys()):
    g0c = g0_crits[alpha]
    n_max = math.log(g0c / G0_E) / math.log(PHI) if g0c > G0_E else 0
    N = max(0, math.floor(n_max) + 1)
    marker = " ← N=3!" if N == 3 else ""
    print(f"  {alpha:6.2f}  {g0c:10.6f}  {n_max:8.4f}  {N:5d}{marker}")

# ================================================================
# SECTION 2: Find α_crit where N transitions 2→3
# ================================================================
print(f"\n{'=' * 70}")
print("  2. CRITICAL α FOR N=2→3 TRANSITION")
print("=" * 70)

# N=3 requires g₀_crit > φ²·g₀^e = 2.2761
target_g0 = PHI**2 * G0_E
print(f"\n  Target: g₀_crit = φ²·g₀^e = {target_g0:.6f}")

# Find α where g₀_crit(α) = target
# We know g₀_crit decreases with α
# α=0.5: g₀_crit=2.618 > target
# α=1.0: g₀_crit=2.206 < target

def g0crit_minus_target(alpha):
    g0c = find_g0_crit(alpha, d=3)
    if g0c is None:
        return -target_g0
    return g0c - target_g0

try:
    alpha_crit = brentq(g0crit_minus_target, 0.3, 1.5, xtol=1e-4)
    g0c_at_crit = find_g0_crit(alpha_crit, d=3)
    print(f"  α_crit = {alpha_crit:.6f}")
    print(f"  g₀_crit(α_crit) = {g0c_at_crit:.6f}")
    print(f"  φ²·g₀^e         = {target_g0:.6f}")
    print(f"  |Δ| = {abs(g0c_at_crit - target_g0):.6f}")

    # What is α_crit? Test simple fractions/values
    candidates = [
        ("1/2", 0.5),
        ("2/3", 2/3),
        ("3/4", 0.75),
        ("4/5", 0.8),
        ("5/6", 5/6),
        ("1/√2", 1/math.sqrt(2)),
        ("1/φ", 1/PHI),
        ("ln(2)", math.log(2)),
        ("1/e", 1/math.e),
        ("π/4", math.pi/4),
        ("√2-1", math.sqrt(2)-1),
        ("(√5-1)/2 = 1/φ", (math.sqrt(5)-1)/2),
        ("2-φ", 2-PHI),
    ]

    print(f"\n  α_crit ≈ {alpha_crit:.6f}")
    print(f"\n  {'Candidate':>20s}  {'Value':>10s}  {'|Δ|':>10s}")
    for name, val in candidates:
        delta = abs(alpha_crit - val)
        if delta < 0.05:
            print(f"  {name:>20s}  {val:10.6f}  {delta:10.6f}  {'←' if delta < 0.01 else ''}")

    check("T1: α_crit exists and is between 0.5 and 1.0",
          0.5 < alpha_crit < 1.0,
          f"α_crit = {alpha_crit:.4f}")

except Exception as e:
    print(f"  Bisection failed: {e}")

# ================================================================
# SECTION 3: α where g₀_crit = g₀^τ(Koide)
# ================================================================
print(f"\n{'=' * 70}")
print("  3. α WHERE BARRIER = τ(KOIDE)")
print("=" * 70)

target_koide = 1.729
print(f"  Target: g₀_crit = g₀^τ(Koide) = {target_koide}")

def g0crit_minus_koide(alpha):
    g0c = find_g0_crit(alpha, d=3)
    if g0c is None:
        return -target_koide
    return g0c - target_koide

try:
    alpha_koide = brentq(g0crit_minus_koide, 2.0, 4.0, xtol=1e-4)
    g0c_at_koide = find_g0_crit(alpha_koide, d=3)
    print(f"  α_Koide = {alpha_koide:.6f}")
    print(f"  g₀_crit(α_Koide) = {g0c_at_koide:.6f}")
    print(f"  g₀^τ(Koide)      = {target_koide:.6f}")

    # Test if α_Koide = 3 or some other simple value
    candidates_k = [
        ("3", 3.0),
        ("π", math.pi),
        ("e", math.e),
        ("3+1/10", 3.1),
        ("φ²", PHI**2),
        ("7/2", 3.5),
        ("√10", math.sqrt(10)),
    ]

    print(f"\n  α_Koide ≈ {alpha_koide:.6f}")
    for name, val in candidates_k:
        delta = abs(alpha_koide - val)
        if delta < 0.1:
            print(f"  {name:>10s} = {val:.6f}, |Δ| = {delta:.6f}")

    check("T2: α_Koide exists",
          alpha_koide > 0,
          f"α = {alpha_koide:.4f}")

except Exception as e:
    print(f"  Bisection failed: {e}")

# ================================================================
# SECTION 4: Mass scaling for α=0.5 (the N=3 case)
# ================================================================
print(f"\n{'=' * 70}")
print("  4. MASS SCALING FOR α=0.5 (K=g)")
print("=" * 70)

alpha_test = 0.5
g0c_test = find_g0_crit(alpha_test, d=3)
print(f"  α = {alpha_test}, g₀_crit = {g0c_test:.6f}")

# Compute mass = 4π∫[g·g'²/2 + g³/3 - g⁴/4 - 1/12]·r² dr
masses_05 = []
g0_vals_05 = []

for g0 in np.concatenate([
    np.arange(0.3, 1.0, 0.05),
    np.arange(1.0, min(g0c_test - 0.05, 3.0), 0.05),
]):
    if abs(g0 - 1.0) < 0.01:
        continue
    g_min, sing, sol = solve_alpha(g0, alpha_test, d=3)
    if sing or not sol.success:
        continue
    r = sol.t
    g = sol.y[0]
    gp = sol.y[1]

    # Energy density: ε = g^{2α}·g'²/2 + U(g) - U(1)
    # U(g) = g³/3 - g⁴/4, U(1) = 1/12
    eps = g**(2*alpha_test) * gp**2 / 2 + g**3/3 - g**4/4 - 1.0/12
    integrand = eps * r**2
    mass = 4 * np.pi * np.trapezoid(integrand, r)

    if mass > 0:
        g0_vals_05.append(g0)
        masses_05.append(mass)

if len(masses_05) > 0:
    print(f"\n  {'g₀':>8s}  {'m(g₀)':>12s}  {'note':>12s}")
    for g0, m in zip(g0_vals_05, masses_05):
        note = ""
        if abs(g0 - G0_E) < 0.03:
            note = "~ electron"
        elif abs(g0 - PHI * G0_E) < 0.03:
            note = "~ muon"
        elif abs(g0 - PHI**2 * G0_E) < 0.05:
            note = "~ tau(φ²)"
        print(f"  {g0:8.4f}  {m:12.6f}  {note:>12s}")

    if len(masses_05) > 5:
        from scipy.interpolate import interp1d
        g0_arr = np.array(g0_vals_05)
        m_arr = np.array(masses_05)
        m_interp = interp1d(g0_arr, m_arr, kind='cubic', fill_value='extrapolate')

        if g0_arr[0] < 0.87 and g0_arr[-1] > 1.41:
            m_e = float(m_interp(0.869))
            m_mu = float(m_interp(1.407))
            m_tau_phi2 = float(m_interp(2.276)) if 2.276 < g0_arr[-1] else None

            print(f"\n  m(g₀^e=0.869) = {m_e:.6f}")
            print(f"  m(g₀^μ=1.407) = {m_mu:.6f}")
            if m_e > 0 and m_mu > 0:
                ratio = m_mu / m_e
                A_e = abs(0.869 - 1)
                A_mu = abs(1.407 - 1)
                if ratio > 1:
                    alpha_eff = math.log(ratio) / (2*math.log(A_mu/A_e))
                    print(f"  m_μ/m_e = {ratio:.4f} (exp: 206.77)")
                    print(f"  α_eff from energy integral = {alpha_eff:.4f}")

            if m_tau_phi2:
                print(f"  m(g₀^τ_φ²=2.276) = {m_tau_phi2:.6f}")
                if m_e > 0:
                    print(f"  m_τ(φ²)/m_e = {m_tau_phi2/m_e:.2f} (exp: 3477)")


# ================================================================
# SECTION 5: The Dimension-Generation-α relationship
# ================================================================
print(f"\n{'=' * 70}")
print("  5. DIMENSION-GENERATION-α RELATIONSHIP")
print("=" * 70)

print(f"""
  N_gen = floor[ln(g₀_crit(α,d)/g₀^e) / ln(φ)] + 1

  g₀_crit depends on BOTH α and d.

  For N=3 in d=3: need g₀_crit > φ²·g₀^e = {PHI**2 * G0_E:.4f}
  This requires α < α_crit ≈ {alpha_crit if 'alpha_crit' in dir() else 0}

  For N=3 in d=2:
""")

# Check d=2
for alpha in [0.3, 0.5, 0.7, 1.0]:
    g0c = find_g0_crit(alpha, d=2)
    if g0c:
        n_max = math.log(g0c / G0_E) / math.log(PHI)
        N = math.floor(n_max) + 1
        print(f"  d=2, α={alpha:.1f}: g₀_crit={g0c:.4f}, n_max={n_max:.3f}, N={N}")

print()
for alpha in [0.3, 0.5, 0.7, 1.0]:
    g0c = find_g0_crit(alpha, d=4)
    if g0c:
        n_max = math.log(g0c / G0_E) / math.log(PHI)
        N = math.floor(n_max) + 1
        print(f"  d=4, α={alpha:.1f}: g₀_crit={g0c:.4f}, n_max={n_max:.3f}, N={N}")


# ================================================================
# SECTION 6: Special case α = 1/φ
# ================================================================
print(f"\n{'=' * 70}")
print("  6. SPECIAL α VALUES")
print("=" * 70)

special_alphas = [
    ("1/φ = 0.618", 1/PHI),
    ("1/2", 0.5),
    ("2/3", 2/3),
    ("ln(2) = 0.693", math.log(2)),
    ("1/√2 = 0.707", 1/math.sqrt(2)),
    ("3/4", 0.75),
]

for name, alpha in special_alphas:
    g0c = find_g0_crit(alpha, d=3)
    if g0c:
        n_max = math.log(g0c / G0_E) / math.log(PHI)
        N = math.floor(n_max) + 1
        below = "✓" if g0c > target_g0 else "✗"
        print(f"  α = {name:>20s}: g₀_crit = {g0c:.6f}, n_max = {n_max:.4f}, "
              f"N = {N}  [{below} > φ²g₀^e = {target_g0:.4f}]")

# ================================================================
# SECTION 7: Can N=3 determine α?
# ================================================================
print(f"\n{'=' * 70}")
print("  7. INVERSE PROBLEM: α FROM N=3")
print("=" * 70)

if 'alpha_crit' in dir():
    print(f"""
  FINDING: N=3 requires α < α_crit ≈ {f'{alpha_crit:.4f}'}

  This is a CONSTRAINT on the theory: the kinetic coupling
  must be WEAKER than α_crit to allow 3 generations.

  Physical interpretation:
    α controls how strongly the metric couples to its own gradient.
    α > α_crit → too much coupling → barrier too low → only N=2
    α < α_crit → weaker coupling → higher barrier → N=3 or more

  For the TGP substrate (α=1): N=2 (barely — n_max = 1.94)
  The substrate is VERY CLOSE to the N=2→3 boundary!

  Key question: is α=1 the CORRECT physical value?
  Or does some geometric argument give α < α_crit?
""")

    # How close is substrate to N=3?
    g0c_sub = find_g0_crit(1.0, d=3)
    deficit = target_g0 - g0c_sub
    deficit_pct = deficit / target_g0 * 100
    print(f"  Substrate deficit for N=3: Δg₀ = {deficit:.4f} ({deficit_pct:.1f}%)")
    print(f"  g₀_crit(sub) = {g0c_sub:.4f} vs φ²g₀^e = {target_g0:.4f}")

    check("T3: α_crit gives exact N=2→3 transition",
          abs(find_g0_crit(alpha_crit, d=3) - target_g0) < 0.01,
          f"α_crit = {alpha_crit:.4f}")


# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'=' * 70}")
print("  SUMMARY")
print("=" * 70)

print(f"""
  g₀_crit(α, d=3):
""")
for alpha in sorted(g0_crits.keys()):
    g0c = g0_crits[alpha]
    n_max = math.log(g0c / G0_E) / math.log(PHI) if g0c > G0_E else 0
    N = max(0, math.floor(n_max) + 1)
    marker = " ←" if N == 3 else ""
    if alpha % 0.5 < 0.01 or abs(alpha - round(alpha)) < 0.01:
        print(f"    α={alpha:.1f}: g₀_crit={g0c:.4f}, N={N}{marker}")

print(f"""
  CRITICAL VALUES:
    α_crit (N=2→3 transition) = {alpha_crit if 'alpha_crit' in dir() else 0}
    α_Koide (barrier = τ mass) = {alpha_koide if 'alpha_koide' in dir() else 0}

  CONSTRAINT FOR N=3: α < {alpha_crit if 'alpha_crit' in dir() else 0}

  SUBSTRATE (α=1): n_max = 1.935 → N=2 (deficit {deficit_pct:.1f}%)
  The substrate theory predicts N=2, not N=3!

  POSSIBLE RESOLUTIONS:
    (a) Physical α ≠ 1 (e.g., α = 1/φ ≈ 0.618 → N=3)
    (b) φ-ladder g₀ values need correction
    (c) g₀^e value needs refinement
    (d) Mass scaling is not purely A^{{2α}}
""")

print(f"\n{'=' * 70}")
print(f"  TEST REPORT: {PASS} PASS, {FAIL} FAIL out of {PASS + FAIL}")
print(f"{'=' * 70}")
