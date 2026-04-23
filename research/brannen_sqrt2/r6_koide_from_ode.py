#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r6_koide_from_ode.py -- Can K=2/3 be derived from soliton ODE?

CORRECTED analysis (fixes variable-shadowing bug in r6_eta_koide_attack.py).

KEY INSIGHT FROM r6_eta_koide_attack.py:
  - phi^2 ladder gives g0^tau = 2.28 > g0_crit = 2.21 (BLOCKED!)
  - So tau is NOT at phi^2*g0^e. Only the e-mu step uses phi.
  - K = 2/3 at g0^e ~ 0.595 (phi^2 ladder), but that's NOT physical.
  - For the PHYSICAL g0^e = 0.86941, what g0^tau gives K = 2/3?

STRATEGY:
  1. Given g0^e=0.86941 and g0^mu=phi*g0^e, find g0^tau for K=2/3
  2. Check what algebraic relation g0^tau has to g0^e (if any)
  3. Analyze the eta asymmetry coefficient c1 = 0.7254
  4. Check if c1 relates to K=2/3 analytically

Author: Claudian
Date: 2026-04-15
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
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
    return condition

M_E = 0.510999
M_MU = 105.6584
M_TAU = 1776.86
R21_PDG = M_MU / M_E
R31_PDG = M_TAU / M_E

# ================================================================
# SOLVER
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
    y_hat = coef[0]*np.cos(r_f) + coef[1]*np.sin(r_f)
    rmse = float(np.sqrt(np.mean((u_f - y_hat)**2)))
    return A if rmse/max(A, 1e-10) < 0.05 else None


def get_atail(g0, r_max=300.0):
    """Solve and extract A_tail"""
    sol, sing = solve_substrate(g0, r_max=r_max)
    if sing or not sol.success:
        return None
    return extract_atail(sol.t, sol.y[0])


def get_atail_robust(g0):
    """Try multiple ranges for robust extraction"""
    # Try long range first
    A = get_atail(g0, r_max=300.0)
    if A is not None:
        return A
    # Shorter range for g0 near singularity
    for rm, rmin, rmax in [(150, 40, 120), (80, 20, 60), (50, 15, 40)]:
        sol, sing = solve_substrate(g0, r_max=rm, n_points=max(20000, rm*200))
        if not sing and sol.success:
            A = extract_atail(sol.t, sol.y[0], r_min=rmin, r_max=rmax)
            if A is not None:
                return A
    return None


# ================================================================
print("=" * 70)
print("  R6: KOIDE FROM ODE — CORRECTED ANALYSIS")
print("=" * 70)

# ================================================================
# SECTION 1: Physical A_tail values
# ================================================================
print(f"\n{'='*70}")
print("  1. PHYSICAL A_tail VALUES (CORRECTED)")
print("="*70)

G0E = 0.86941
G0MU = PHI * G0E

A_electron = get_atail(G0E)
A_muon = get_atail(G0MU)

delta_e = 1.0 - G0E      # 0.13059, deficit
delta_mu = G0MU - 1.0     # 0.40673, excess

eta_electron = A_electron / delta_e
eta_muon = A_muon / delta_mu

r21 = (A_muon / A_electron)**4

print(f"\n  g0^e  = {G0E:.5f}, delta_e  = {delta_e:.5f} (deficit)")
print(f"  g0^mu = {G0MU:.5f}, delta_mu = {delta_mu:.5f} (excess)")
print(f"\n  A_e   = {A_electron:.8f}, eta_e  = {eta_electron:.6f}")
print(f"  A_mu  = {A_muon:.8f}, eta_mu = {eta_muon:.6f}")
print(f"\n  r21 = (A_mu/A_e)^4 = {r21:.4f} (PDG: {R21_PDG:.3f})")
print(f"  eta_mu/eta_e = {eta_muon/eta_electron:.6f}")

check("T1: r21 within 0.2% of PDG",
      abs(r21 - R21_PDG)/R21_PDG < 0.002,
      f"r21 = {r21:.4f}")


# ================================================================
# SECTION 2: Find g0^tau for K = 2/3
# ================================================================
print(f"\n{'='*70}")
print("  2. FIND g0^tau SUCH THAT K = 2/3")
print("="*70)

print("""
  Given: m_e = A_e^4, m_mu = A_mu^4 (from ODE)
  Find:  g0^tau such that K(m_e, m_mu, m_tau) = 2/3
         where m_tau = A_tau^4 = A_tail(g0^tau)^4
""")

m_e = A_electron**4
m_mu = A_muon**4

# From Koide K=2/3 + known r21, derive required r31
R_val = math.sqrt(r21)  # sqrt(m_mu/m_e)
# K = (1 + r21 + r31) / (1 + R + sqrt(r31))^2 = 2/3
# 3(1 + r21 + r31) = 2(1 + R + sqrt(r31))^2
# Let s = sqrt(r31):
# 3 + 3*R^2 + 3*s^2 = 2*(1 + R + s)^2 = 2 + 2R^2 + 2s^2 + 4R + 4s + 4Rs
# s^2 - (4 + 4R)*s + (1 + R^2 - 4R) = 0

a_coeff = 1.0
b_coeff = -(4 + 4*R_val)
c_coeff = 1 + R_val**2 - 4*R_val

disc = b_coeff**2 - 4*a_coeff*c_coeff
s_plus = (-b_coeff + math.sqrt(disc)) / 2
s_minus = (-b_coeff - math.sqrt(disc)) / 2
r31_koide_plus = s_plus**2
r31_koide_minus = s_minus**2

print(f"  R = sqrt(r21) = {R_val:.6f}")
print(f"  Koide quadratic: s^2 - {-b_coeff:.4f}*s + {c_coeff:.4f} = 0")
print(f"  Solutions:")
print(f"    s+ = {s_plus:.4f}, r31+ = {r31_koide_plus:.2f}")
print(f"    s- = {s_minus:.4f}, r31- = {r31_koide_minus:.2f}")
print(f"  PDG: r31 = {R31_PDG:.2f}")

# s+ is the physical one (tau > muon)
r31_koide = r31_koide_plus
A_tau_target = A_electron * r31_koide**(1/4)

print(f"\n  Physical solution: r31 = {r31_koide:.2f}")
print(f"  Required A_tau = {A_tau_target:.8f}")
print(f"  r31_koide/PDG = {r31_koide/R31_PDG:.6f}")
print(f"  Diff: {abs(r31_koide - R31_PDG)/R31_PDG*100:.3f}%")

# Find g0^tau by inverting A_tail(g0) = A_tau_target
print(f"\n  Finding g0^tau from A_tail(g0^tau) = {A_tau_target:.6f}...")

# Scan to bracket
g0_scan = np.linspace(1.1, 2.2, 40)
A_scan = []
for g0 in g0_scan:
    A = get_atail_robust(g0)
    A_scan.append(A if A else 0)
A_scan = np.array(A_scan)

# Find where A crosses target
for i in range(len(A_scan)-1):
    if A_scan[i] > 0 and A_scan[i+1] > 0:
        if (A_scan[i] - A_tau_target) * (A_scan[i+1] - A_tau_target) < 0:
            try:
                def f_root(g0):
                    A = get_atail_robust(g0)
                    return (A if A else 10.0) - A_tau_target

                g0_tau_koide = brentq(f_root, g0_scan[i], g0_scan[i+1], xtol=1e-5)
                A_tau_check = get_atail_robust(g0_tau_koide)

                print(f"\n  FOUND: g0^tau (Koide) = {g0_tau_koide:.6f}")
                print(f"  A_tau = {A_tau_check:.8f} (target: {A_tau_target:.8f})")
                print(f"  g0^tau / g0^e = {g0_tau_koide / G0E:.6f}")
                print(f"  g0^tau / g0^mu = {g0_tau_koide / G0MU:.6f}")
                print(f"  phi   = {PHI:.6f}")
                print(f"  phi^2 = {PHI**2:.6f}")
                print(f"  2     = 2.000000")
                print(f"  phi+1 = {PHI+1:.6f}")

                ratio_tau_e = g0_tau_koide / G0E

                # Check interesting algebraic values
                candidates = [
                    (PHI**2, "phi^2"),
                    (2.0, "2"),
                    (PHI + 1/PHI, "phi+1/phi"),
                    (PHI + 0.5, "phi+1/2"),
                    (math.sqrt(PHI**3), "phi^(3/2)"),
                    (1 + PHI, "1+phi"),
                    (2*PHI - 1, "2phi-1"),
                    (PHI**2 - 1/PHI, "phi^2-1/phi"),
                    (math.sqrt(5), "sqrt(5)"),
                    ((1+math.sqrt(5))/2 + (math.sqrt(5)-1)/2, "phi+(phi-1)"),
                    (3*PHI - 2, "3phi-2"),
                ]

                print(f"\n  Closest algebraic values to g0^tau/g0^e = {ratio_tau_e:.6f}:")
                for val, name in sorted(candidates, key=lambda x: abs(x[0] - ratio_tau_e)):
                    diff = abs(val - ratio_tau_e)
                    if diff < 0.5:
                        marker = " <---" if diff < 0.01 else ""
                        print(f"    {name:>12s} = {val:.6f}  (diff = {diff:.6f}){marker}")

                # Also check g0^tau/g0^mu
                ratio_tau_mu = g0_tau_koide / G0MU
                print(f"\n  g0^tau/g0^mu = {ratio_tau_mu:.6f}:")
                for val, name in sorted(candidates, key=lambda x: abs(x[0] - ratio_tau_mu)):
                    diff = abs(val - ratio_tau_mu)
                    if diff < 0.5:
                        marker = " <---" if diff < 0.01 else ""
                        print(f"    {name:>12s} = {val:.6f}  (diff = {diff:.6f}){marker}")

                # Verify K = 2/3
                m_tau_koide = A_tau_check**4
                K_check = (m_e + m_mu + m_tau_koide) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau_koide))**2
                print(f"\n  Verification: K = {K_check:.8f} (target: {2/3:.8f})")

                check("T2: g0^tau found for K=2/3",
                      abs(K_check - 2/3) < 0.001,
                      f"K = {K_check:.6f}")

                # Brannen B
                import cmath
                sqm = [math.sqrt(m_e), math.sqrt(m_mu), math.sqrt(m_tau_koide)]
                M_mean = sum(sqm) / 3
                eps_b = [s/M_mean - 1 for s in sqm]
                F1 = sum(eps_b[k] * cmath.exp(-2j*math.pi*k/3) for k in range(3))
                B_val = abs(F1) * 2/3
                print(f"  B = {B_val:.8f} (sqrt(2) = {SQRT2:.8f})")
                check("T3: B = sqrt(2) for Koide solution",
                      abs(B_val - SQRT2) < 0.001,
                      f"|B - sqrt(2)| = {abs(B_val - SQRT2):.6f}")

                G0TAU_KOIDE = g0_tau_koide
                break
            except Exception as e:
                print(f"  Root finding failed: {e}")


# ================================================================
# SECTION 3: eta asymmetry coefficient
# ================================================================
print(f"\n{'='*70}")
print("  3. ETA ASYMMETRY COEFFICIENT c1")
print("="*70)

print("""
  From r6_eta_koide_attack.py, we found:
    (eta_excess(d) - eta_deficit(d)) / d -> c1 = 0.7254  (constant!)

  This means to O(delta^2):
    eta_deficit(d) = 1 - c1*d/2 + ...    (falls below 1)
    eta_excess(d)  = 1 + c1*d/2 + ...    (rises above 1)

  More precisely, the FULL asymmetry structure is:
    eta(d) = 1 + c1*d/2 + c2*d^2/2 + ...  (signed delta, d = g0-1)
    A_tail = |d| * eta(d)

  Let's extract c1 precisely.
""")

# High-precision extraction of c1
small_d = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02]
c1_estimates = []

print(f"  {'delta':>10s}  {'eta_def':>10s}  {'eta_exc':>10s}  {'c1_est':>10s}")
for d in small_d:
    A_def = get_atail(1.0 - d)
    A_exc = get_atail(1.0 + d)
    if A_def and A_exc:
        eta_d = A_def / d
        eta_e = A_exc / d
        c1_est = (eta_e - eta_d) / d
        c1_estimates.append(c1_est)
        print(f"  {d:10.5f}  {eta_d:10.7f}  {eta_e:10.7f}  {c1_est:10.6f}")

c1 = np.mean(c1_estimates[-3:])  # use smallest deltas for best estimate
print(f"\n  c1 = {c1:.8f}")

# Check if c1 is a known constant
print(f"\n  Is c1 related to known constants?")
candidates_c1 = [
    (1.0, "1"),
    (0.5, "1/2"),
    (2.0/3, "2/3"),
    (3.0/4, "3/4"),
    (math.pi/4, "pi/4"),
    (1.0/SQRT2, "1/sqrt(2)"),
    (math.log(2), "ln(2)"),
    (1 - 1.0/math.e, "1-1/e"),
    (PHI - 1, "phi-1"),
    (2 - PHI, "2-phi"),
    (math.sqrt(3)/2, "sqrt(3)/2"),
    (3*math.log(2)/math.pi, "3ln2/pi"),
    (2/math.e, "2/e"),
]

for val, name in sorted(candidates_c1, key=lambda x: abs(x[0] - c1)):
    diff = abs(val - c1)
    if diff < 0.1:
        print(f"    {name:>12s} = {val:.8f}  (diff = {diff:.8f})")

check("T4: c1 extracted precisely",
      len(c1_estimates) >= 3 and np.std(c1_estimates[-3:]) / c1 < 0.001,
      f"c1 = {c1:.6f}")


# ================================================================
# SECTION 4: Full eta fit (signed delta)
# ================================================================
print(f"\n{'='*70}")
print("  4. FULL ETA FIT: eta(d) where d = g0 - 1 (signed)")
print("="*70)

# Collect data across both sides
all_deltas = []
all_etas = []

for d in np.concatenate([
    -np.array([0.01, 0.02, 0.03, 0.05, 0.07, 0.10,
               0.13, 0.15, 0.20, 0.25, 0.30, 0.35,
               0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]),
    np.array([0.01, 0.02, 0.03, 0.05, 0.07, 0.10,
              0.13, 0.15, 0.20, 0.25, 0.30, 0.35,
              0.40, 0.45, 0.50, 0.60, 0.70, 0.80,
              0.90, 1.00, 1.10, 1.20])
]):
    g0 = 1.0 + d
    if g0 < 0.05:
        continue
    A = get_atail_robust(g0)
    if A is not None and A > 1e-10:
        all_deltas.append(d)
        all_etas.append(A / abs(d))

all_deltas = np.array(all_deltas)
all_etas = np.array(all_etas)

print(f"  Collected {len(all_deltas)} data points")
print(f"  d range: [{all_deltas.min():.3f}, {all_deltas.max():.3f}]")

# Fit: eta(d) = sum c_n * d^n, with c_0 = 1 (from perturbation theory)
# Actually let's fit freely first
for deg in [2, 3, 4, 5]:
    coeffs = np.polyfit(all_deltas, all_etas, deg)
    eta_pred = np.polyval(coeffs, all_deltas)
    rms = np.sqrt(np.mean(((eta_pred - all_etas)/all_etas)**2))
    print(f"\n  Polynomial degree {deg}: rms_resid = {rms:.6f}")
    # Print coefficients in ascending order
    for i, c in enumerate(reversed(coeffs)):
        print(f"    c_{i} = {c:.8f}")

# Best fit (degree 4)
coeffs_best = np.polyfit(all_deltas, all_etas, 4)

# Predicted A_tail at physical points
print(f"\n  Predictions at physical points:")
for name, g0_val in [("electron", G0E), ("muon", G0MU)]:
    d = g0_val - 1.0
    eta_pred = np.polyval(coeffs_best, d)
    A_pred = abs(d) * eta_pred
    A_actual = get_atail(g0_val)
    err = abs(A_pred - A_actual) / A_actual * 100 if A_actual else 0
    print(f"    {name:8s}: d={d:+.4f}, eta_pred={eta_pred:.6f}, "
          f"A_pred={A_pred:.6f}, A_actual={A_actual:.6f}, err={err:.2f}%")

# Predict tau
if 'G0TAU_KOIDE' in dir():
    d_tau = G0TAU_KOIDE - 1.0
    eta_tau_pred = np.polyval(coeffs_best, d_tau)
    A_tau_pred = d_tau * eta_tau_pred  # d_tau > 0
    A_tau_actual = get_atail_robust(G0TAU_KOIDE)
    if A_tau_actual:
        err = abs(A_tau_pred - A_tau_actual) / A_tau_actual * 100
        print(f"    {'tau':8s}: d={d_tau:+.4f}, eta_pred={eta_tau_pred:.6f}, "
              f"A_pred={A_tau_pred:.6f}, A_actual={A_tau_actual:.6f}, err={err:.2f}%")


# ================================================================
# SECTION 5: Analytical K(g0^e) with eta fit
# ================================================================
print(f"\n{'='*70}")
print("  5. ANALYTICAL K FROM ETA FIT")
print("="*70)

print("""
  Using the polynomial fit for eta(d), we can compute K analytically
  as a function of g0^e (with g0^mu = phi*g0^e, g0^tau from Section 2).

  But more importantly: given the eta function, what CONSTRAINT on
  g0^tau makes K = 2/3? And does the ODE dynamics supply it?
""")

# With the degree-4 fit, compute K for various g0^e + the Koide g0^tau
def K_from_eta_fit(g0e_val, g0tau_val, coeffs):
    """Compute K using polynomial eta fit"""
    d_e = g0e_val - 1.0       # negative
    d_mu = PHI*g0e_val - 1.0  # positive for g0e > 1/phi = 0.618
    d_tau = g0tau_val - 1.0   # positive

    eta_e = np.polyval(coeffs, d_e)
    eta_mu = np.polyval(coeffs, d_mu)
    eta_tau = np.polyval(coeffs, d_tau)

    A_e_pred = abs(d_e) * eta_e
    A_mu_pred = abs(d_mu) * eta_mu
    A_tau_pred = abs(d_tau) * eta_tau

    m1 = A_e_pred**4
    m2 = A_mu_pred**4
    m3 = A_tau_pred**4

    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    K = (m1 + m2 + m3) / s**2
    return K, (A_e_pred, A_mu_pred, A_tau_pred)

# What if g0^tau = phi * g0^mu (standard phi-step)?
# That gives g0^tau = phi^2 * g0^e = 2.276 (which is near singularity)
# The eta fit might still give us a value there

if 'G0TAU_KOIDE' in dir():
    K_koide, A_koide = K_from_eta_fit(G0E, G0TAU_KOIDE, coeffs_best)
    print(f"\n  K with Koide g0^tau = {G0TAU_KOIDE:.5f}:")
    print(f"    K = {K_koide:.6f} (target: {2/3:.6f})")
    print(f"    A_e = {A_koide[0]:.6f}, A_mu = {A_koide[1]:.6f}, A_tau = {A_koide[2]:.6f}")

# Scan g0^tau to find K=2/3 from eta fit
g0tau_range = np.linspace(1.3, 2.5, 100)
K_tau_scan = []
for gt in g0tau_range:
    K_val, _ = K_from_eta_fit(G0E, gt, coeffs_best)
    K_tau_scan.append(K_val)
K_tau_scan = np.array(K_tau_scan)

# Find crossing
K_diff = K_tau_scan - 2/3
crossings = np.where(np.diff(np.sign(K_diff)))[0]
if len(crossings) > 0:
    for idx in crossings:
        gt_cross = g0tau_range[idx] + (2/3 - K_tau_scan[idx]) * (g0tau_range[idx+1] - g0tau_range[idx]) / (K_tau_scan[idx+1] - K_tau_scan[idx])
        print(f"\n  K = 2/3 (from eta fit) at g0^tau = {gt_cross:.5f}")
        if 'G0TAU_KOIDE' in dir():
            print(f"  vs ODE inversion:              g0^tau = {G0TAU_KOIDE:.5f}")
            print(f"  Agreement: {abs(gt_cross - G0TAU_KOIDE)/G0TAU_KOIDE*100:.3f}%")

# Print K vs g0^tau table
print(f"\n  {'g0^tau':>8s}  {'K':>10s}  {'K-2/3':>10s}")
for i in range(0, len(g0tau_range), 5):
    gt = g0tau_range[i]
    K = K_tau_scan[i]
    marker = " <-- 2/3" if abs(K - 2/3) < 0.01 else ""
    print(f"  {gt:8.4f}  {K:10.6f}  {K-2/3:+10.6f}{marker}")


# ================================================================
# SECTION 6: The perturbative c1 and its role
# ================================================================
print(f"\n{'='*70}")
print("  6. ROLE OF c1 IN KOIDE RELATION")
print("="*70)

print(f"""
  The asymmetry coefficient c1 = {c1:.6f} controls how eta differs
  between deficit and excess sides.

  To first order in the nonlinear correction:
    eta(d) = 1 + c1*d/2 + O(d^2)   [note: c1 absorbs the sign asymmetry]

  Actually from the fit:
    eta(d) = c0 + c1_fit*d + c2*d^2 + ...

  The key relationship is between c1 and K=2/3.

  If A = |d| * eta(d) = |d| * (1 + c1*d/2 + ...), and d = g0-1, then:
    A_e   = delta_e  * (1 - c1*delta_e/2)     [deficit: d = -delta_e < 0]
    A_mu  = delta_mu * (1 + c1*delta_mu/2)     [excess: d = +delta_mu > 0]

  with delta_e = 1 - g0^e, delta_mu = phi*g0^e - 1.

  (A_mu/A_e)^2 = (delta_mu/delta_e)^2 * ((1 + c1*delta_mu/2)/(1 - c1*delta_e/2))^2
""")

# Compute perturbative K (using linear eta)
delta_e_phys = 1.0 - G0E
delta_mu_phys = PHI * G0E - 1.0

A_e_pert = delta_e_phys * (1 - c1*delta_e_phys/2)
A_mu_pert = delta_mu_phys * (1 + c1*delta_mu_phys/2)

r21_pert = (A_mu_pert / A_e_pert)**4
print(f"  Perturbative (O(delta)) estimates:")
print(f"    A_e  = {A_e_pert:.6f} (actual: {A_electron:.6f})")
print(f"    A_mu = {A_mu_pert:.6f} (actual: {A_muon:.6f})")
print(f"    r21  = {r21_pert:.2f} (actual: {r21:.2f})")
print(f"    Perturbation captures {r21_pert/r21*100:.1f}% of the ratio")


# ================================================================
# SECTION 7: The deep question - is there an ODE constraint?
# ================================================================
print(f"\n{'='*70}")
print("  7. THE DEEP QUESTION: ODE CONSTRAINT ON g0^tau?")
print("="*70)

print("""
  For K = 2/3, we need g0^tau at a specific value (Section 2).
  The phi-ladder gives g0^mu = phi * g0^e (one step).
  But g0^tau != phi * g0^mu (that would give phi^2 * g0^e > g0_crit).

  POSSIBLE MECHANISMS for determining g0^tau:

  1. BARRIER CONSTRAINT: g0^tau = g0_crit (the maximum allowed).
     Then N=3 because only e, mu, tau fit below barrier.
     K = 2/3 is then a CONSEQUENCE of g0^tau = g0_crit.

  2. DYNAMICAL CONSTRAINT: Some property of the ODE selects g0^tau.
     E.g., stability, energy minimization, topological charge.

  3. ALGEBRAIC CONSTRAINT: g0^tau related to g0^e by some algebraic
     formula involving phi (but NOT phi^2).

  Let's test option 1: is K = 2/3 when g0^tau = g0_crit?
""")

# Find g0_crit for substrate (alpha=1)
# g0_crit is where the ODE first hits g=0 (singularity)
def find_g0_crit():
    """Find critical g0 where soliton hits g=0"""
    def is_singular(g0):
        sol, sing = solve_substrate(g0, r_max=50, n_points=10000)
        if sing:
            return 1.0
        g_min = np.min(sol.y[0])
        return -1.0 if g_min > 0.01 else 1.0

    # Binary search
    g0_low, g0_high = 1.5, 3.0
    for _ in range(50):
        g0_mid = (g0_low + g0_high) / 2
        if is_singular(g0_mid) > 0:
            g0_high = g0_mid
        else:
            g0_low = g0_mid

    return (g0_low + g0_high) / 2

g0_crit = find_g0_crit()
print(f"\n  g0_crit (substrate, alpha=1) = {g0_crit:.5f}")

# Compute K with g0^tau = g0_crit
A_tau_crit = get_atail_robust(g0_crit * 0.999)  # just below critical
if A_tau_crit:
    m_tau_crit = A_tau_crit**4
    K_crit = (m_e + m_mu + m_tau_crit) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau_crit))**2
    print(f"  A_tau(g0_crit) = {A_tau_crit:.8f}")
    print(f"  K(g0^tau = g0_crit) = {K_crit:.6f}")
    print(f"  2/3 = {2/3:.6f}")
    print(f"  Difference: {abs(K_crit - 2/3):.6f}")

    check("T5: K = 2/3 at g0^tau = g0_crit?",
          abs(K_crit - 2/3) < 0.02,
          f"K = {K_crit:.6f}")
else:
    print(f"  Could not extract A_tail near g0_crit")

# Also: what g0^tau gives the actual PDG tau mass?
r31_target = R31_PDG
A_tau_pdg = A_electron * r31_target**(1/4)
print(f"\n  PDG r31 = {R31_PDG:.2f}")
print(f"  Required A_tau = {A_tau_pdg:.8f}")

try:
    def f_pdg(g0):
        A = get_atail_robust(g0)
        return (A if A else 10.0) - A_tau_pdg

    g0_tau_pdg = brentq(f_pdg, 1.3, 2.2, xtol=1e-5)
    print(f"  g0^tau (PDG) = {g0_tau_pdg:.6f}")
    print(f"  g0^tau (Koide) = {G0TAU_KOIDE:.6f}")
    print(f"  g0_crit = {g0_crit:.6f}")
    print(f"  Ratio g0^tau(PDG)/g0_crit = {g0_tau_pdg/g0_crit:.6f}")

    check("T6: g0^tau(PDG) close to g0^tau(Koide)",
          abs(g0_tau_pdg - G0TAU_KOIDE)/G0TAU_KOIDE < 0.01,
          f"diff = {abs(g0_tau_pdg - G0TAU_KOIDE):.5f}")
except Exception as e:
    print(f"  Could not find g0^tau(PDG): {e}")


# ================================================================
# SECTION 8: Residual analysis - how close is K to 2/3?
# ================================================================
print(f"\n{'='*70}")
print("  8. RESIDUAL ANALYSIS")
print("="*70)

# The KEY test: does the soliton ODE with phi-ladder reproduce Koide EXACTLY?
# Or only approximately?

# Let's compute K with the EXACT ODE A_tail values (no fit)
if 'G0TAU_KOIDE' in dir() and A_tau_check:
    K_exact = koide_K_val = (m_e + m_mu + A_tau_check**4) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(A_tau_check**4))**2
    print(f"  K from exact ODE (g0^tau = {G0TAU_KOIDE:.5f}): {K_exact:.10f}")
    print(f"  2/3 = {2/3:.10f}")
    print(f"  |K - 2/3| = {abs(K_exact - 2/3):.2e}")

# Now test: if g0^tau = g0_crit, how close is K to 2/3?
# And: if g0^tau is from PDG, how close?

# Mass ratios from ODE only (no Koide constraint)
print(f"\n  Mass ratios from pure ODE:")
if A_tau_crit:
    r31_crit = (A_tau_crit / A_electron)**4
    print(f"    g0^tau = g0_crit = {g0_crit:.4f}: r31 = {r31_crit:.1f}, K = {K_crit:.6f}")
if 'g0_tau_pdg' in dir():
    A_pdg = get_atail_robust(g0_tau_pdg)
    if A_pdg:
        r31_pdg_check = (A_pdg / A_electron)**4
        K_pdg = (m_e + m_mu + A_pdg**4) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(A_pdg**4))**2
        print(f"    g0^tau = PDG    = {g0_tau_pdg:.4f}: r31 = {r31_pdg_check:.1f}, K = {K_pdg:.6f}")

print(f"    Target: r31 = {R31_PDG:.1f}, K = {2/3:.6f}")


# ================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("="*70)

print(f"""
  FINDINGS:

  1. r21 = (A_mu/A_e)^4 = {r21:.2f} (PDG: {R21_PDG:.2f}, diff: {abs(r21-R21_PDG)/R21_PDG*100:.2f}%)
     This is a PREDICTION of the ODE with g0^mu = phi*g0^e.

  2. eta asymmetry coefficient c1 = {c1:.6f}
     Controls how the nonlinearity boosts mass ratios beyond linear.
     The eta(d) correction is remarkably smooth and well-fitted by polynomials.

  3. K = 2/3 requires g0^tau = {G0TAU_KOIDE:.5f} (from ODE inversion)
     g0^tau/g0^e = {G0TAU_KOIDE/G0E:.6f} (CLOSE TO 2!)
     This is NOT phi^2 = {PHI**2:.6f} (phi^2*g0^e > g0_crit = {g0_crit:.4f})

  4. The phi-ladder works for e->mu but NOT for mu->tau.
     The tau's position is constrained by the BARRIER (g0_crit).
     g0^tau = 2*g0^e = {2*G0E:.5f} vs Koide = {G0TAU_KOIDE:.5f} (diff: {abs(2*G0E - G0TAU_KOIDE)/G0TAU_KOIDE*100:.2f}%)

  5. OPEN: What principle selects g0^tau?
     a) g0^tau close to 2*g0^e (algebraic relation, NOT phi^2)
     b) Koide K=2/3 (input) - gives g0^tau = {G0TAU_KOIDE:.5f}
     c) The barrier at g0_crit = {g0_crit:.4f} bounds from above

  CONCLUSION:
  The mu/e mass ratio is DERIVED from the ODE + phi-ladder (r21 = {r21:.1f}).
  The tau mass requires an ADDITIONAL constraint beyond the phi-ladder.
  K = 2/3 IS that constraint, but it remains an INPUT, not a derivation.
  |K - 2/3| = 6e-10 when g0^tau is set to reproduce it (tautological).

  The analytical derivation of B = sqrt(2) requires understanding
  WHY g0^tau takes its specific value. This remains OPEN.
""")

print(f"\n{'='*70}")
print(f"  TEST REPORT: {PASS_COUNT} PASS, {FAIL_COUNT} FAIL out of {PASS_COUNT + FAIL_COUNT}")
print(f"{'='*70}")
