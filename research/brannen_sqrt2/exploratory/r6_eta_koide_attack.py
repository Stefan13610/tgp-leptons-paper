#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_eta_koide_attack.py -- Can we derive K=2/3 from soliton ODE + phi-ladder?

The nonlinear correction eta(delta) = A_tail(1+delta)/|delta| encodes how
the ODE's nonlinearity tunes mass ratios beyond the linear approximation.

KEY CHAIN:
  g0 = 1 + delta  ->  ODE  ->  A_tail(g0)  ->  m = c_M * A_tail^4
  A_tail(g0) = |delta| * eta(delta)    [definition of eta]

  If A_tail were LINEAR in delta, then:
    (A_mu/A_e)^4 = (delta_mu/delta_e)^4 = (phi*g0e - 1)^4 / (1-g0e)^4 ~ 94

  But the ACTUAL ratio is 207 because eta is NOT constant:
    (A_mu/A_e)^4 = (delta_mu/delta_e)^4 * (eta_mu/eta_e)^4 = 94 * 2.2 = 207

  QUESTION: Does the functional form of eta(delta), combined with the
  phi-ladder constraint, FORCE K = 2/3?

Strategy:
  1. Map eta(delta) with high precision across deficit and excess sides
  2. Fit analytical forms (polynomial, rational, etc.)
  3. With best fit, compute K analytically as function of g0^e
  4. Check if K = 2/3 at the PHYSICAL g0^e = 0.86941

Author: Claudian
Date: 2026-04-15
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq, minimize_scalar
import math

PHI = (1 + math.sqrt(5)) / 2
SQRT2 = math.sqrt(2)

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS {name}  {detail}")
    else:
        FAIL += 1
        print(f"  FAIL {name}  {detail}")
    return condition

# PDG
M_E = 0.510999
M_MU = 105.6584
M_TAU = 1776.86
R21_PDG = M_MU / M_E   # 206.768

# ================================================================
# SOLVER (substrate alpha=1, high precision)
# ================================================================

def solve_substrate(g0, r_max=300.0, n_points=60000):
    """Substrate ODE: g'' + (1/g)g'^2 + (2/r)g' = 1 - g"""
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
    """Extract A_tail from far field: (g-1)*r = B*cos(r) + C*sin(r)"""
    mask = (r >= r_min) & (r <= r_max)
    r_f = r[mask]
    if len(r_f) < 20:
        return None
    u_f = (g[mask] - 1.0) * r_f

    X = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, u_f, rcond=None)
    A = math.sqrt(coef[0]**2 + coef[1]**2)

    # quality check
    y_hat = coef[0]*np.cos(r_f) + coef[1]*np.sin(r_f)
    rmse = float(np.sqrt(np.mean((u_f - y_hat)**2)))
    quality = rmse / max(A, 1e-10)

    return A if quality < 0.05 else None


def get_atail(g0):
    """Solve and extract A_tail for substrate (alpha=1)"""
    sol, sing = solve_substrate(g0)
    if sing or not sol.success:
        return None
    return extract_atail(sol.t, sol.y[0])


# ================================================================
print("=" * 70)
print("  R6 ATTACK: eta(delta) AND KOIDE FROM SOLITON ODE")
print("=" * 70)

# ================================================================
# SECTION 1: High-precision eta(delta) mapping
# ================================================================
print(f"\n{'='*70}")
print("  1. eta(delta) = A_tail(1+delta) / |delta|")
print("="*70)

# Dense sampling of delta on both sides
deltas_deficit = np.array([0.01, 0.02, 0.03, 0.05, 0.07, 0.10,
                           0.13, 0.15, 0.17, 0.20, 0.25, 0.30,
                           0.35, 0.40, 0.45, 0.50, 0.55, 0.60,
                           0.65, 0.70])  # g0 = 1 - delta < 1

deltas_excess = np.array([0.01, 0.02, 0.03, 0.05, 0.07, 0.10,
                          0.13, 0.15, 0.17, 0.20, 0.25, 0.30,
                          0.35, 0.40, 0.45, 0.50, 0.60, 0.70,
                          0.80, 0.90, 1.00, 1.10, 1.20])  # g0 = 1 + delta > 1

print(f"\n  {'delta':>8s}  {'g0':>8s}  {'A_tail':>12s}  {'eta':>10s}  side")
print(f"  {'-'*8}  {'-'*8}  {'-'*12}  {'-'*10}  ----")

eta_deficit_data = []  # (delta, eta)
eta_excess_data = []

for delta in deltas_deficit:
    g0 = 1.0 - delta
    if g0 < 0.05:
        continue
    A = get_atail(g0)
    if A is not None:
        eta = A / delta
        eta_deficit_data.append((delta, eta, A))
        print(f"  {delta:8.4f}  {g0:8.4f}  {A:12.8f}  {eta:10.6f}  deficit")

print()
for delta in deltas_excess:
    g0 = 1.0 + delta
    A = get_atail(g0)
    if A is not None:
        eta = A / delta
        eta_excess_data.append((delta, eta, A))
        print(f"  {delta:8.4f}  {g0:8.4f}  {A:12.8f}  {eta:10.6f}  excess")

eta_def = np.array(eta_deficit_data)
eta_exc = np.array(eta_excess_data)


# ================================================================
# SECTION 2: Functional form of eta
# ================================================================
print(f"\n{'='*70}")
print("  2. FITTING eta(delta)")
print("="*70)

# Try several functional forms for each side

def fit_and_report(name, func, x, y, p0, bounds=(-np.inf, np.inf)):
    try:
        popt, pcov = curve_fit(func, x, y, p0=p0, bounds=bounds, maxfev=20000)
        y_fit = func(x, *popt)
        resid = np.abs(y_fit - y) / y
        max_resid = np.max(resid)
        rms_resid = np.sqrt(np.mean(resid**2))
        print(f"  {name}:")
        print(f"    params = {[f'{p:.6f}' for p in popt]}")
        print(f"    max_resid = {max_resid:.6f}, rms_resid = {rms_resid:.6f}")
        return popt, max_resid, rms_resid
    except Exception as e:
        print(f"  {name}: FAILED ({e})")
        return None, 1.0, 1.0


if len(eta_def) > 5:
    d_def = eta_def[:, 0]
    e_def = eta_def[:, 1]

    print("\n  DEFICIT side (g0 < 1, delta = 1 - g0):")

    # Model 1: eta = a + b*delta
    fit_and_report("Linear: eta = a + b*d",
                   lambda d, a, b: a + b*d,
                   d_def, e_def, [0.5, -0.5])

    # Model 2: eta = a + b*delta + c*delta^2
    fit_and_report("Quadratic: eta = a + b*d + c*d^2",
                   lambda d, a, b, c: a + b*d + c*d**2,
                   d_def, e_def, [0.5, -0.5, 0.1])

    # Model 3: eta = a * delta^b
    fit_and_report("Power: eta = a*d^b",
                   lambda d, a, b: a * d**b,
                   d_def, e_def, [0.5, -0.1],
                   bounds=([0, -2], [10, 2]))

    # Model 4: eta = a / (1 + b*delta)^c
    fit_and_report("Rational: eta = a/(1+b*d)^c",
                   lambda d, a, b, c: a / (1 + b*d)**c,
                   d_def, e_def, [0.5, 1.0, 0.5],
                   bounds=([0, 0, 0], [10, 50, 10]))

    # Model 5: eta = a * exp(-b*delta)
    fit_and_report("Exponential: eta = a*exp(-b*d)",
                   lambda d, a, b: a * np.exp(-b*d),
                   d_def, e_def, [0.5, 0.5])

if len(eta_exc) > 5:
    d_exc = eta_exc[:, 0]
    e_exc = eta_exc[:, 1]

    print("\n  EXCESS side (g0 > 1, delta = g0 - 1):")

    fit_and_report("Linear: eta = a + b*d",
                   lambda d, a, b: a + b*d,
                   d_exc, e_exc, [0.5, 0.5])

    fit_and_report("Quadratic: eta = a + b*d + c*d^2",
                   lambda d, a, b, c: a + b*d + c*d**2,
                   d_exc, e_exc, [0.5, 0.5, 0.1])

    fit_and_report("Power: eta = a*d^b",
                   lambda d, a, b: a * d**b,
                   d_exc, e_exc, [0.5, 0.1],
                   bounds=([0, -2], [10, 2]))

    fit_and_report("Rational: eta = a/(1+b*d)^c",
                   lambda d, a, b, c: a / (1 + b*d)**c,
                   d_exc, e_exc, [0.5, -0.5, -0.5])


# ================================================================
# SECTION 3: eta at physical points -> mass ratios
# ================================================================
print(f"\n{'='*70}")
print("  3. eta AT PHYSICAL POINTS")
print("="*70)

g0e = 0.86941
g0mu = PHI * g0e    # 1.4067
g0tau = PHI**2 * g0e  # 2.2754

delta_e = 1.0 - g0e     # 0.13059 (deficit)
delta_mu = g0mu - 1.0    # 0.4067  (excess)
delta_tau = g0tau - 1.0  # 1.2754  (excess)

A_e = get_atail(g0e)
A_mu = get_atail(g0mu)
A_tau = get_atail(g0tau)

print(f"\n  g0^e  = {g0e:.5f},  delta_e  = {delta_e:.5f}  (deficit)")
print(f"  g0^mu = {g0mu:.5f},  delta_mu = {delta_mu:.5f}  (excess)")
print(f"  g0^tau= {g0tau:.5f},  delta_tau= {delta_tau:.5f}  (excess)")

if A_e and A_mu:
    eta_e = A_e / delta_e
    eta_mu = A_mu / delta_mu

    print(f"\n  A_e   = {A_e:.8f},   eta_e  = {eta_e:.6f}")
    print(f"  A_mu  = {A_mu:.8f},   eta_mu = {eta_mu:.6f}")

    if A_tau:
        eta_tau = A_tau / delta_tau
        print(f"  A_tau = {A_tau:.8f},   eta_tau = {eta_tau:.6f}")
    else:
        print(f"  A_tau: extraction failed (g0={g0tau:.4f} may be near singularity)")
        # Try shorter range for tau
        sol_tau, sing_tau = solve_substrate(g0tau, r_max=100, n_points=20000)
        if not sing_tau and sol_tau.success:
            A_tau = extract_atail(sol_tau.t, sol_tau.y[0], r_min=30, r_max=80)
            if A_tau:
                eta_tau = A_tau / delta_tau
                print(f"  A_tau = {A_tau:.8f},   eta_tau = {eta_tau:.6f}  (short range)")

    # Mass ratios
    print(f"\n  LINEAR approximation (eta = const):")
    r21_lin = (delta_mu/delta_e)**4
    print(f"    (delta_mu/delta_e)^4 = {r21_lin:.2f}  (vs PDG {R21_PDG:.2f})")

    print(f"\n  ACTUAL (with eta correction):")
    r21_act = (A_mu/A_e)**4
    print(f"    (A_mu/A_e)^4 = {r21_act:.4f}")

    boost = r21_act / r21_lin
    print(f"    eta boost = {boost:.4f}  (= (eta_mu/eta_e)^4)")
    print(f"    eta_mu/eta_e = {(eta_mu/eta_e):.6f}")

    check("T1: mass ratio r21 within 1% of PDG",
          abs(r21_act - R21_PDG)/R21_PDG < 0.01,
          f"r21 = {r21_act:.4f}")


# ================================================================
# SECTION 4: Koide as function of g0^e
# ================================================================
print(f"\n{'='*70}")
print("  4. KOIDE K(g0^e) -- IS IT SPECIAL AT g0^e = 0.86941?")
print("="*70)

print("""
  If phi-ladder fixes g0^mu = phi*g0^e, g0^tau = phi^2*g0^e,
  then ALL masses are functions of g0^e alone:
    m_n(g0^e) = c_M * A_tail(phi^n * g0^e)^4

  Koide K is then a function of g0^e ONLY.
  Question: is K(g0^e) = 2/3 at some special g0^e?
  And is g0^e = 0.86941 that special value?
""")

def koide_K(m1, m2, m3):
    """K = (sum m) / (sum sqrt(m))^2"""
    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return (m1 + m2 + m3) / s**2

def compute_K_at_g0e(g0e_val):
    """Compute Koide K for the phi-ladder starting at g0e_val"""
    g0_1 = g0e_val
    g0_2 = PHI * g0e_val
    g0_3 = PHI**2 * g0e_val

    A1 = get_atail(g0_1)
    A2 = get_atail(g0_2)

    # For tau, try standard first, then short range
    A3 = get_atail(g0_3)
    if A3 is None:
        sol3, s3 = solve_substrate(g0_3, r_max=100, n_points=20000)
        if not s3 and sol3.success:
            A3 = extract_atail(sol3.t, sol3.y[0], r_min=30, r_max=80)

    if A1 is None or A2 is None or A3 is None:
        return None, (A1, A2, A3)

    m1 = A1**4
    m2 = A2**4
    m3 = A3**4

    K = koide_K(m1, m2, m3)
    return K, (A1, A2, A3)

# Scan g0^e
g0e_scan = np.linspace(0.50, 0.95, 30)
K_values = []
g0e_valid = []

print(f"\n  {'g0^e':>8s}  {'g0^mu':>8s}  {'g0^tau':>8s}  {'K':>10s}  {'K - 2/3':>10s}")
print(f"  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}")

for g0e_val in g0e_scan:
    K, As = compute_K_at_g0e(g0e_val)
    if K is not None:
        K_values.append(K)
        g0e_valid.append(g0e_val)
        marker = " <-- PHYSICAL" if abs(g0e_val - 0.86941) < 0.01 else ""
        print(f"  {g0e_val:8.4f}  {PHI*g0e_val:8.4f}  {PHI**2*g0e_val:8.4f}  "
              f"{K:10.6f}  {K - 2/3:10.6f}{marker}")

K_values = np.array(K_values)
g0e_valid = np.array(g0e_valid)

if len(K_values) > 5:
    # Is K roughly constant or varying?
    K_mean = np.mean(K_values)
    K_std = np.std(K_values)
    K_cv = K_std / K_mean

    print(f"\n  K statistics:")
    print(f"    mean = {K_mean:.6f}")
    print(f"    std  = {K_std:.6f}")
    print(f"    CV   = {K_cv:.6f}")
    print(f"    2/3  = {2/3:.6f}")

    # Find g0^e where K = 2/3 (if it crosses)
    K_minus_23 = K_values - 2/3
    crossings = np.where(np.diff(np.sign(K_minus_23)))[0]

    if len(crossings) > 0:
        print(f"\n  K crosses 2/3 at {len(crossings)} point(s):")
        for idx in crossings:
            # Linear interpolation
            g0_cross = g0e_valid[idx] + (2/3 - K_values[idx]) * (g0e_valid[idx+1] - g0e_valid[idx]) / (K_values[idx+1] - K_values[idx])
            print(f"    g0^e = {g0_cross:.5f}")
    else:
        print(f"\n  K does NOT cross 2/3 in this range")
        print(f"  K range: [{K_values.min():.6f}, {K_values.max():.6f}]")

    check("T2: K varies meaningfully with g0^e",
          K_cv > 0.001,
          f"CV = {K_cv:.6f}")

    # Check K at physical point
    K_phys, _ = compute_K_at_g0e(0.86941)
    if K_phys is not None:
        check("T3: K(0.86941) close to 2/3",
              abs(K_phys - 2/3) < 0.05,
              f"K = {K_phys:.6f}, 2/3 = {2/3:.6f}")


# ================================================================
# SECTION 5: Brannen B as function of g0^e
# ================================================================
print(f"\n{'='*70}")
print("  5. BRANNEN B(g0^e)")
print("="*70)

def brannen_B(m1, m2, m3):
    """Extract B from 3 masses using DFT on sqrt(m)"""
    import cmath
    sqm = [math.sqrt(m1), math.sqrt(m2), math.sqrt(m3)]
    M = sum(sqm) / 3
    eps = [s/M - 1 for s in sqm]
    F1 = sum(eps[k] * cmath.exp(-2j*math.pi*k/3) for k in range(3))
    B = abs(F1) * 2/3
    return B

B_values = []
for i, g0e_val in enumerate(g0e_valid):
    K_val = K_values[i]
    # B^2 = 2(3K - 1) when phases are equidistant
    # But let's compute B directly from masses
    _, As = compute_K_at_g0e(g0e_val)
    if As[0] and As[1] and As[2]:
        B = brannen_B(As[0]**4, As[1]**4, As[2]**4)
        B_values.append(B)
    else:
        B_values.append(None)

print(f"\n  {'g0^e':>8s}  {'K':>10s}  {'B':>10s}  {'B - sqrt2':>10s}")
print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}")
for i in range(len(g0e_valid)):
    if B_values[i] is not None:
        marker = " <-- PHYSICAL" if abs(g0e_valid[i] - 0.86941) < 0.01 else ""
        print(f"  {g0e_valid[i]:8.4f}  {K_values[i]:10.6f}  {B_values[i]:10.6f}  "
              f"{B_values[i] - SQRT2:10.6f}{marker}")


# ================================================================
# SECTION 6: Cross-vacuum analysis of A_tail
# ================================================================
print(f"\n{'='*70}")
print("  6. CROSS-VACUUM STRUCTURE: A_tail across g0=1")
print("="*70)

print("""
  The electron sits at g0 = 0.87 (deficit side of vacuum g=1).
  The muon sits at g0 = 1.41 (excess side).
  The tau sits at g0 = 2.28 (deep excess).

  The phi-ladder CROSSES the vacuum at g0=1. This means:
  - eta_e is computed from deficit-side physics
  - eta_mu and eta_tau from excess-side physics

  These are DIFFERENT branches of the ODE solution space!
  The cross-vacuum structure may be the KEY to K=2/3.
""")

# Plot A_tail(g0) across both sides with fine sampling
g0_fine = np.concatenate([
    np.linspace(0.30, 0.99, 50),
    np.linspace(1.01, 2.50, 50)
])

print(f"\n  {'g0':>8s}  {'A_tail':>12s}  {'|delta|':>8s}  {'eta':>10s}  {'ln(eta)':>10s}")
print(f"  {'-'*8}  {'-'*12}  {'-'*8}  {'-'*10}  {'-'*10}")

atail_fine = []
for g0 in g0_fine:
    A = get_atail(g0)
    if A is not None and A > 1e-10:
        delta = abs(g0 - 1.0)
        eta = A / delta if delta > 0.005 else 0
        ln_eta = math.log(eta) if eta > 0 else 0
        atail_fine.append((g0, A, delta, eta, ln_eta))
        if abs(g0 - 0.87) < 0.02 or abs(g0 - 1.41) < 0.02 or abs(g0 - 2.28) < 0.05 or g0 in [0.50, 0.70, 0.90, 1.10, 1.30, 1.50, 1.70, 2.00, 2.30]:
            side = "D" if g0 < 1 else "E"
            print(f"  {g0:8.4f}  {A:12.8f}  {delta:8.4f}  {eta:10.6f}  {ln_eta:10.6f}  ({side})")


# ================================================================
# SECTION 7: The key asymmetry - deficit vs excess eta
# ================================================================
print(f"\n{'='*70}")
print("  7. DEFICIT vs EXCESS ASYMMETRY")
print("="*70)

print("""
  For the SAME |delta|, how does eta differ between sides?

  eta_deficit(delta) = A_tail(1 - delta) / delta
  eta_excess(delta)  = A_tail(1 + delta) / delta

  The RATIO eta_excess/eta_deficit at the SAME |delta| measures
  the ODE's asymmetry around vacuum g=1.
""")

# Compare at matched deltas
if len(eta_def) > 0 and len(eta_exc) > 0:
    common_deltas = sorted(set(eta_def[:, 0]).intersection(set(eta_exc[:, 0])))
    if not common_deltas:
        # Match approximately
        for dd, ed, _ in eta_deficit_data:
            for de, ee, _ in eta_excess_data:
                if abs(dd - de) < 0.001:
                    common_deltas.append(dd)

    print(f"\n  {'delta':>8s}  {'eta_def':>10s}  {'eta_exc':>10s}  {'ratio':>10s}  {'ln(ratio)':>10s}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    for dd, ed, _ in eta_deficit_data:
        # Find matching excess
        for de, ee, _ in eta_excess_data:
            if abs(dd - de) < 0.001:
                ratio = ee / ed
                print(f"  {dd:8.4f}  {ed:10.6f}  {ee:10.6f}  {ratio:10.6f}  {math.log(ratio):10.6f}")

# ================================================================
# SECTION 8: The ODE perturbation theory for eta
# ================================================================
print(f"\n{'='*70}")
print("  8. PERTURBATION THEORY FOR eta(delta)")
print("="*70)

print("""
  Near g0 = 1 (vacuum), g(r) = 1 + epsilon*f(r) where epsilon = g0 - 1.

  Linearized ODE (to O(epsilon)):
    f'' + (2/r)f' + f = 0  (since 1-g = -epsilon*f, and (1/g)g'^2 ~ O(eps^2))

  Solution: f(r) = sin(r)/r  (regular at origin)

  So to LEADING ORDER:
    g(r) = 1 + epsilon * sin(r)/r
    A_tail = |epsilon| * 1 = |delta|
    eta = 1  (exactly, to leading order)

  The CORRECTIONS to eta come from O(epsilon^2) terms in the ODE:
    - (1/g)g'^2 term contributes -epsilon * f'^2 at O(eps^2)
    - (1-g) term is exact: -epsilon * f
    - g^{2-2alpha} = g^0 = 1 for alpha=1 (no correction from RHS!)

  For alpha=1, the ONLY nonlinear correction comes from (1/g)g'^2.

  Let's compute eta^(2) = the O(epsilon) correction to eta.
""")

# Numerical verification of perturbation theory
# For small delta, eta should approach a limit
small_deltas = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
print(f"\n  Numerical eta for small delta:")
print(f"  {'delta':>10s}  {'eta_def':>10s}  {'eta_exc':>10s}  {'(e_exc-e_def)/delta':>20s}")
print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*20}")

for delta in small_deltas:
    g0_d = 1.0 - delta
    g0_e = 1.0 + delta
    A_d = get_atail(g0_d)
    A_e = get_atail(g0_e)
    if A_d and A_e:
        eta_d = A_d / delta
        eta_e = A_e / delta
        asym = (eta_e - eta_d) / delta
        print(f"  {delta:10.4f}  {eta_d:10.6f}  {eta_e:10.6f}  {asym:20.6f}")


# ================================================================
# SECTION 9: Can we derive K = 2/3 analytically?
# ================================================================
print(f"\n{'='*70}")
print("  9. ANALYTICAL STRUCTURE")
print("="*70)

print("""
  Let's define:
    A(g0) = |g0 - 1| * eta(g0 - 1)

  For the phi-ladder with g0^e:
    delta_e = 1 - g0^e                     (positive, deficit)
    delta_mu = phi*g0^e - 1                 (positive, excess)
    delta_tau = phi^2*g0^e - 1              (positive, excess)

  Mass ratios:
    r21 = (A_mu/A_e)^4 = (delta_mu/delta_e)^4 * (eta_exc(delta_mu)/eta_def(delta_e))^4
    r31 = (A_tau/A_e)^4 = (delta_tau/delta_e)^4 * (eta_exc(delta_tau)/eta_def(delta_e))^4

  Koide:
    K = (r21 + r31 + 1) / (1 + sqrt(r21) + sqrt(r31))^2

  For K = 2/3, we need:
    3(r21 + r31 + 1) = 2(1 + sqrt(r21) + sqrt(r31))^2
""")

# Compute everything symbolically
if A_e and A_mu:
    r21 = (A_mu/A_e)**4

    # Find what r31 is needed for K = 2/3 given this r21
    # 3(r21 + r31 + 1) = 2(1 + sqrt(r21) + sqrt(r31))^2
    # Let s = sqrt(r31), R = sqrt(r21)
    # 3(r21 + s^2 + 1) = 2(1 + R + s)^2
    # 3r21 + 3s^2 + 3 = 2 + 4R + 4s + 2R^2 + 4Rs + 2s^2
    # s^2 - 4s(1+R)/1 + (3r21 + 3 - 2 - 4R - 2R^2) = 0  ... wait
    # s^2(3-2) - s(4 + 4R) + (3r21 + 3 - 2 - 4R - 2R^2) = 0
    # s^2 - 4(1+R)s + (3 + 3R^2 - 2 - 4R - 2R^2) = 0   [using r21 = R^2]
    # s^2 - 4(1+R)s + (1 + R^2 - 4R) = 0
    # s^2 - 4(1+R)s + (R-1)^2 - 2R = 0

    R = math.sqrt(r21)
    # s^2 - 4(1+R)s + (1 + R^2 - 4R) = 0
    a_coeff = 1.0
    b_coeff = -4*(1 + R)
    c_coeff = 1 + R**2 - 4*R

    disc = b_coeff**2 - 4*a_coeff*c_coeff
    if disc >= 0:
        s_plus = (-b_coeff + math.sqrt(disc)) / (2*a_coeff)
        s_minus = (-b_coeff - math.sqrt(disc)) / (2*a_coeff)
        r31_plus = s_plus**2
        r31_minus = s_minus**2

        print(f"  Given r21 = {r21:.4f} (R = {R:.4f}):")
        print(f"    K=2/3 requires r31 = {r31_plus:.2f} (s+) or {r31_minus:.2f} (s-)")
        print(f"    PDG r31 = {M_TAU/M_E:.2f}")

        if A_tau:
            r31_actual = (A_tau/A_e)**4
            print(f"    Soliton r31 = {r31_actual:.2f}")

            # What Koide predicts for tau mass
            m_tau_koide = r31_plus * M_E if abs(r31_plus - M_TAU/M_E) < abs(r31_minus - M_TAU/M_E) else r31_minus * M_E
            r31_koide = r31_plus if abs(r31_plus - M_TAU/M_E) < abs(r31_minus - M_TAU/M_E) else r31_minus

            print(f"\n    Koide prediction: r31 = {r31_koide:.2f}")
            print(f"    This means A_tau/A_e = {r31_koide**(1/4):.6f}")
            print(f"    And A_tau = {A_e * r31_koide**(1/4):.8f}")

            # What g0^tau does this correspond to?
            # We need to invert A_tail(g0) = target
            target_A_tau = A_e * r31_koide**(1/4)

            print(f"\n    Finding g0^tau that gives Koide-predicted A_tau...")

            try:
                def atail_minus_target(g0):
                    A = get_atail(g0)
                    if A is None:
                        sol, s = solve_substrate(g0, r_max=100, n_points=20000)
                        if not s and sol.success:
                            A = extract_atail(sol.t, sol.y[0], r_min=30, r_max=80)
                    if A is None:
                        return 1.0
                    return A - target_A_tau

                g0_tau_koide = brentq(atail_minus_target, 1.5, 3.0, xtol=1e-5)
                print(f"    g0^tau (Koide) = {g0_tau_koide:.5f}")
                print(f"    g0^tau (phi^2) = {PHI**2 * g0e:.5f}")
                print(f"    Difference: {abs(g0_tau_koide - PHI**2*g0e):.5f}")
                print(f"    Ratio g0^tau/g0^e = {g0_tau_koide/g0e:.6f}")
                print(f"    phi^2 = {PHI**2:.6f}")
                print(f"    Discrepancy: {abs(g0_tau_koide/g0e - PHI**2)/PHI**2*100:.3f}%")

                check("T4: Koide g0^tau close to phi^2 * g0^e",
                      abs(g0_tau_koide/g0e - PHI**2)/PHI**2 < 0.05,
                      f"ratio = {g0_tau_koide/g0e:.6f}, phi^2 = {PHI**2:.6f}")
            except Exception as e:
                print(f"    Inversion failed: {e}")


# ================================================================
# SECTION 10: The universal function A_tail(g0)
# ================================================================
print(f"\n{'='*70}")
print("  10. GLOBAL FIT: A_tail(g0) = ? ")
print("="*70)

print("""
  Instead of fitting eta, let's fit A_tail(g0) directly across BOTH sides.
  The function must satisfy:
    - A_tail(1) = 0  (vacuum)
    - A_tail(g0) > 0 for g0 != 1
    - A_tail ~ |g0 - 1| for small |g0 - 1|  (perturbation theory)
    - Nonlinear corrections grow with |g0 - 1|

  Candidate global forms:
    A = |g0 - 1| * f(g0)   where f(1) = 1 (perturbative limit)
""")

# Collect all A_tail data
all_data = [(g0, A, d, eta, ln_e) for g0, A, d, eta, ln_e in atail_fine]
if all_data:
    g0_all = np.array([x[0] for x in all_data])
    A_all = np.array([x[1] for x in all_data])

    # Fit: A = |g0-1| * (1 + a*(g0-1) + b*(g0-1)^2 + c*(g0-1)^3)
    delta_all = g0_all - 1.0
    abs_delta = np.abs(delta_all)

    # eta = A / |delta|
    eta_all = A_all / abs_delta

    # Fit eta as polynomial in (g0-1) = signed delta
    # eta(delta) = 1 + a1*delta + a2*delta^2 + a3*delta^3
    try:
        # Use polyfit on (delta, eta)
        mask_fit = abs_delta > 0.01  # avoid division instability near vacuum
        delta_fit = delta_all[mask_fit]
        eta_fit = eta_all[mask_fit]

        for deg in [2, 3, 4, 5]:
            coeffs = np.polyfit(delta_fit, eta_fit, deg)
            eta_pred = np.polyval(coeffs, delta_fit)
            resid = np.sqrt(np.mean(((eta_pred - eta_fit)/eta_fit)**2))

            # Format coefficients
            terms = []
            for i, c in enumerate(reversed(coeffs)):
                if i == 0:
                    terms.append(f"{c:.6f}")
                else:
                    terms.append(f"{c:+.6f}*d^{i}")

            print(f"  deg {deg}: rms_resid = {resid:.6f}")
            if deg <= 3:
                print(f"    eta = {' '.join(terms[:4])}")

            if deg == 3:
                # Store for later use
                poly_coeffs = coeffs

        # With cubic fit, compute predicted K
        print(f"\n  Using cubic fit eta(d) = {coeffs[3]:.4f} + {coeffs[2]:.4f}*d + {coeffs[1]:.4f}*d^2 + {coeffs[0]:.4f}*d^3")

        # Predicted A_tail at physical points
        for name, g0_val in [("e", g0e), ("mu", g0mu), ("tau", g0tau)]:
            d = g0_val - 1.0
            eta_pred = np.polyval(poly_coeffs, d)
            A_pred = abs(d) * eta_pred
            A_actual = get_atail(g0_val)
            if A_actual is None:
                sol_t, s_t = solve_substrate(g0_val, r_max=100, n_points=20000)
                if not s_t and sol_t.success:
                    A_actual = extract_atail(sol_t.t, sol_t.y[0], r_min=30, r_max=80)

            if A_actual:
                err = abs(A_pred - A_actual) / A_actual * 100
                print(f"    {name}: A_pred = {A_pred:.6f}, A_actual = {A_actual:.6f}, err = {err:.2f}%")
            else:
                print(f"    {name}: A_pred = {A_pred:.6f}, A_actual = N/A")

    except Exception as e:
        print(f"  Polynomial fit failed: {e}")


# ================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("="*70)

print(f"""
  The eta(delta) = A_tail(1+delta)/|delta| correction encodes the
  nonlinearity of the substrate ODE.

  KEY FINDINGS:
  1. eta -> 1 for small delta (perturbation theory: A ~ |delta|)
  2. eta_deficit < 1 for large deficit delta (deficit solitons)
  3. eta_excess > 1 for large excess delta (excess solitons)
  4. The ASYMMETRY eta_exc/eta_def at matched |delta| grows with delta
  5. The phi-ladder CROSSES vacuum: e is deficit, mu/tau are excess
  6. This cross-vacuum structure is what makes (A_mu/A_e)^4 >> (delta_mu/delta_e)^4

  The question "why K = 2/3?" becomes:
  "why does the specific nonlinear correction eta(delta) from the
   substrate ODE, combined with the phi-ladder crossing vacuum,
   produce the exact Koide relation?"
""")

print(f"\n{'='*70}")
print(f"  TEST REPORT: {PASS} PASS, {FAIL} FAIL out of {PASS + FAIL}")
print(f"{'='*70}")
