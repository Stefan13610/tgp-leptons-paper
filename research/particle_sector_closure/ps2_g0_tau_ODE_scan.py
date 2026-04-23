#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ps2_g0_tau_ODE_scan.py  --  Program P4, problem #2.

GOAL:  First-principles determination of g_0^tau.  Currently g_0^tau = 1.729
       is FIT to Koide K = 2/3.  Can it be derived instead from:
         (i)  PDG mass ratio m_tau/m_e = 3477.15 via  m_i/m_j = (A_i/A_j)^4, OR
         (ii) A joint constraint (Koide + mass + sum=4) that uniquely picks it?

SETUP:
  Substrate ODE (alpha=1, d=3):
     g'' + (1/g)(g')^2 + (2/r) g' = (1-g) g^{2-2*1} = (1-g) g^0 = 1 - g

  Wait: in brannen_sqrt2 notation the ODE is (see r3_el_check.py,
  r3_conservation_universal.py):
     g'' + (alpha/g)(g')^2 + ((d-1)/r) g' = (1-g) * g^{2-2*alpha}

  For alpha=1, d=3:
     g'' + (1/g)(g')^2 + (2/r) g' = (1-g) * g^0 = 1 - g

  Boundary conditions:
     g(0) = g_0   (central value, parameter)
     g'(0) = 0    (regularity)
     g(r -> inf) -> 1  (asymptotic vacuum)

  For large r: g - 1 = A * sin(r + delta) / r   (linearized tail)

STRATEGY:
  (A) Solve ODE for g_0 in [0.5, 2.2] with 20+ points, extract A_tail(g_0).
  (B) Bisect g_0^tau so that (A_tau/A_e)^4 = m_tau/m_e = 3477.15 (PDG).
  (C) Compute Koide K with this g_0^tau;  compare to 2/3.
  (D) Independently bisect g_0^tau so that K_Koide = 2/3 exactly.
  (E) Compare the two g_0^tau values.  If they agree to <0.01%, Koide is
      DERIVED from ODE + PDG masses.  If they disagree, the residual is
      a prediction (testable vs PDG).
  (F) Test the sum(g_0) = 4 conservation law:  g_0^e + g_0^mu + g_0^tau = 4?

Runtime: ~30 seconds (several ODE solves with RK45).
"""
import sys
import io
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

print("=" * 78)
print("  ps2_g0_tau_ODE_scan.py")
print("=" * 78)
print()
print("  Problem P4 #2:  derive g_0^tau from first principles.")
print()

# ---------------------------------------------------------------------------
# Part A.  Substrate ODE solver
# ---------------------------------------------------------------------------
PHI = (1 + np.sqrt(5)) / 2

# PDG mass ratios (2024)
M_E = 0.51099895069  # MeV
M_MU = 105.6583755
M_TAU = 1776.86
R21_PDG = M_MU / M_E          # 206.7683
R31_PDG = M_TAU / M_E         # 3477.15
K_KOIDE_TARGET = 2.0 / 3.0


def rhs(r, y):
    g, gp = y
    # Avoid singular g=0 (never reached physically in our range)
    if g < 1e-12:
        g = 1e-12
    # ODE:  g'' = (1 - g) - (1/g) (g')^2 - (2/r) g'
    # At r=0 use L'Hopital:  (2/r) g' -> 2 g''(0)  =>  g''(0) = (1-g0)/3 - (g'0)^2/(3g0)
    if r < 1e-8:
        gpp = (1 - g) / 3.0 - 0.0  # g'(0)=0
    else:
        gpp = (1 - g) - (gp * gp) / g - 2.0 * gp / r
    return [gp, gpp]


def solve_soliton(g0, r_max=80.0, rtol=1e-10, atol=1e-12):
    """Solve ODE from r=0 to r=r_max."""
    y0 = [g0, 0.0]
    t_span = (1e-6, r_max)
    # Small shift from r=0 with Taylor expansion
    r_start = 1e-4
    g_start = g0 + (1 - g0) / 6.0 * r_start ** 2  # Taylor: g(r) ~ g0 + ((1-g0)/6)*r^2
    gp_start = (1 - g0) / 3.0 * r_start           # g'(r) ~ ((1-g0)/3)*r
    y0 = [g_start, gp_start]
    sol = solve_ivp(
        rhs, (r_start, r_max), y0,
        method="DOP853", rtol=rtol, atol=atol, dense_output=True,
    )
    return sol


def extract_A_tail(sol, r_fit_min=25.0, r_fit_max=60.0, n_pts=400):
    """Fit g - 1 = A sin(r + delta)/r over [r_fit_min, r_fit_max].

    Uses linear fit of form (g-1)*r = A cos(delta) sin(r) + A sin(delta) cos(r).
    """
    r_arr = np.linspace(r_fit_min, r_fit_max, n_pts)
    y_arr = sol.sol(r_arr)
    g_arr = y_arr[0]
    rhs_vec = (g_arr - 1.0) * r_arr
    # Basis: sin(r), cos(r)
    M = np.column_stack([np.sin(r_arr), np.cos(r_arr)])
    coef, *_ = np.linalg.lstsq(M, rhs_vec, rcond=None)
    a_sin, a_cos = coef
    A = np.sqrt(a_sin ** 2 + a_cos ** 2)
    delta = np.arctan2(a_cos, a_sin)
    # Check quality
    resid = rhs_vec - M @ coef
    rms = float(np.sqrt(np.mean(resid ** 2)))
    return A, delta, rms


# ---------------------------------------------------------------------------
# Part B.  Scan g_0 -> A_tail, establish A_e baseline
# ---------------------------------------------------------------------------
print("=" * 78)
print("  Part A+B.  Scan g_0 -> A_tail")
print("=" * 78)

# g_0^e from Compton matching (canonical value from TGP)
g0_e = 0.86941
sol_e = solve_soliton(g0_e)
A_e, delta_e, rms_e = extract_A_tail(sol_e)
print(f"  g_0^e    = {g0_e:.6f}   A_e = {A_e:.8f}   rms = {rms_e:.2e}")

# g_0^mu = phi * g_0^e (phi-ladder hypothesis)
g0_mu = PHI * g0_e
sol_mu = solve_soliton(g0_mu)
A_mu, delta_mu, rms_mu = extract_A_tail(sol_mu)
ratio_21 = (A_mu / A_e) ** 4
print(f"  g_0^mu   = phi*g_0^e = {g0_mu:.6f}   A_mu = {A_mu:.8f}")
print(f"  (A_mu/A_e)^4 = {ratio_21:.4f}   PDG r_21 = {R21_PDG:.4f}")
print(f"  diff from PDG: {100*(ratio_21 - R21_PDG)/R21_PDG:+.3f}%")

# ---------------------------------------------------------------------------
# Part C.  Bisect g_0^tau so that (A_tau/A_e)^4 = r_31 (PDG)
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part C.  Bisect g_0^tau via PDG mass:   (A_tau/A_e)^4 = r_31")
print("=" * 78)


def objective_mass(g0_tau):
    sol_tau = solve_soliton(g0_tau)
    A_tau, _, _ = extract_A_tail(sol_tau)
    return (A_tau / A_e) ** 4 - R31_PDG


t0 = time.time()
g0_tau_mass = brentq(objective_mass, 1.4, 2.2, xtol=1e-8, rtol=1e-10)
print(f"  g_0^tau (from PDG m_tau/m_e)  = {g0_tau_mass:.8f}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Part D.  Bisect g_0^tau so that Koide K = 2/3 exactly
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part D.  Bisect g_0^tau via Koide K = 2/3")
print("=" * 78)


def koide_K(m1, m2, m3):
    num = m1 + m2 + m3
    denom = (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)) ** 2
    return num / denom


def objective_koide(g0_tau):
    sol_tau = solve_soliton(g0_tau)
    A_tau, _, _ = extract_A_tail(sol_tau)
    # Masses (arbitrary overall scale; ratios matter for K)
    m_e_rel = A_e ** 4
    m_mu_rel = A_mu ** 4
    m_tau_rel = A_tau ** 4
    return koide_K(m_e_rel, m_mu_rel, m_tau_rel) - K_KOIDE_TARGET


t0 = time.time()
g0_tau_koide = brentq(objective_koide, 1.4, 2.2, xtol=1e-8, rtol=1e-10)
print(f"  g_0^tau (from Koide K=2/3)    = {g0_tau_koide:.8f}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Part E.  Comparison - does ODE + PDG masses DERIVE Koide?
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part E.  Comparison:  mass-derived vs Koide-derived g_0^tau")
print("=" * 78)

diff_tau = g0_tau_mass - g0_tau_koide
diff_pct = 100.0 * diff_tau / g0_tau_koide
print(f"  g_0^tau (mass-derived)   = {g0_tau_mass:.8f}")
print(f"  g_0^tau (Koide-derived)  = {g0_tau_koide:.8f}")
print(f"  diff                      = {diff_tau:+.2e}  ({diff_pct:+.3f}%)")
print()

# Evaluate both at each other's point
sol_mt = solve_soliton(g0_tau_mass)
A_tau_m, _, _ = extract_A_tail(sol_mt)
K_at_mass = koide_K(A_e ** 4, A_mu ** 4, A_tau_m ** 4)
m_at_mass = (A_tau_m / A_e) ** 4

sol_kt = solve_soliton(g0_tau_koide)
A_tau_k, _, _ = extract_A_tail(sol_kt)
K_at_koide = koide_K(A_e ** 4, A_mu ** 4, A_tau_k ** 4)
m_at_koide = (A_tau_k / A_e) ** 4

print(f"  AT g_0^tau (mass-derived) = {g0_tau_mass:.6f}:")
print(f"     (A_tau/A_e)^4 = {m_at_mass:.4f}   [target {R31_PDG:.4f}]")
print(f"     Koide K       = {K_at_mass:.6f}   [target {K_KOIDE_TARGET:.6f}]   "
      f"diff = {100*(K_at_mass-K_KOIDE_TARGET):+.3f}%")
print()
print(f"  AT g_0^tau (Koide-derived) = {g0_tau_koide:.6f}:")
print(f"     (A_tau/A_e)^4 = {m_at_koide:.4f}   [target {R31_PDG:.4f}]   "
      f"diff = {100*(m_at_koide-R31_PDG)/R31_PDG:+.3f}%")
print(f"     Koide K       = {K_at_koide:.6f}   [target {K_KOIDE_TARGET:.6f}]")

# ---------------------------------------------------------------------------
# Part F.  Sum conservation:  g_0^e + g_0^mu + g_0^tau = 4?
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part F.  Sum conservation:   g_0^e + g_0^mu + g_0^tau = 4?")
print("=" * 78)

sum_mass = g0_e + g0_mu + g0_tau_mass
sum_koide = g0_e + g0_mu + g0_tau_koide
print(f"  SUM (mass-derived)  = {g0_e:.4f} + {g0_mu:.4f} + {g0_tau_mass:.4f} = "
      f"{sum_mass:.6f}   diff from 4 = {sum_mass - 4:+.5f}")
print(f"  SUM (Koide-derived) = {g0_e:.4f} + {g0_mu:.4f} + {g0_tau_koide:.4f} = "
      f"{sum_koide:.6f}   diff from 4 = {sum_koide - 4:+.5f}")

# Solve: what would g_0^tau need to be for SUM = 4 exactly?
g0_tau_sum = 4.0 - g0_e - g0_mu
sol_st = solve_soliton(g0_tau_sum)
A_tau_s, _, _ = extract_A_tail(sol_st)
m_at_sum = (A_tau_s / A_e) ** 4
K_at_sum = koide_K(A_e ** 4, A_mu ** 4, A_tau_s ** 4)
print()
print(f"  IF SUM = 4 exactly:  g_0^tau = 4 - g_0^e - g_0^mu = {g0_tau_sum:.6f}")
print(f"     (A_tau/A_e)^4 = {m_at_sum:.4f}   [PDG {R31_PDG:.4f}]  "
      f"diff = {100*(m_at_sum-R31_PDG)/R31_PDG:+.3f}%")
print(f"     Koide K       = {K_at_sum:.6f}   [2/3]  "
      f"diff = {100*(K_at_sum-K_KOIDE_TARGET):+.3f}%")

# ---------------------------------------------------------------------------
# Part G.  Algebraic hypotheses on g_0^tau
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part G.  Test algebraic hypotheses on g_0^tau")
print("=" * 78)

hypotheses = [
    ("2 * g_0^e",            2.0 * g0_e),
    ("sqrt(3/2) * g_0^mu",   np.sqrt(1.5) * g0_mu),
    ("phi^2 * g_0^e",        PHI ** 2 * g0_e),
    ("(g_0^e + g_0^mu) * phi / 2",  (g0_e + g0_mu) * PHI / 2),
    ("4 - g_0^e - g_0^mu",   4.0 - g0_e - g0_mu),
    ("sqrt(3) * g_0^e",      np.sqrt(3) * g0_e),
    ("g_0^mu + g_0^e/pi",    g0_mu + g0_e / np.pi),
    ("g_0^mu * (1 + 1/(2 pi))", g0_mu * (1 + 1/(2*np.pi))),
]
print(f"  Reference g_0^tau values:")
print(f"     mass-derived  = {g0_tau_mass:.6f}")
print(f"     Koide-derived = {g0_tau_koide:.6f}")
print()
for name, val in hypotheses:
    d_m = 100 * (val - g0_tau_mass) / g0_tau_mass
    d_k = 100 * (val - g0_tau_koide) / g0_tau_koide
    marker = " <-- CLOSE" if abs(d_k) < 0.2 else ""
    print(f"     {name:35s} = {val:.6f}   diff_mass={d_m:+.3f}%   diff_koide={d_k:+.3f}%{marker}")

# ---------------------------------------------------------------------------
# Part H.  Verdict
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part H.  Verdict")
print("=" * 78)

if abs(diff_pct) < 0.05:
    print("  RESULT:  mass-derived and Koide-derived g_0^tau AGREE to < 0.05%.")
    print("           Koide K=2/3 is DERIVED (not independent of PDG masses).")
    print("           ps2 CLOSURE: the TGP mass formula (m = c K^2, K ~ A^2)")
    print("           together with PDG m_tau/m_e FORCES Koide relation.")
elif abs(diff_pct) < 0.5:
    print(f"  RESULT:  mass-derived and Koide-derived g_0^tau agree to {abs(diff_pct):.2f}%.")
    print("           Consistency good but not perfect.")
    print("           Likely source: fit-window systematic in A_tail extraction.")
    print("           ps2 PARTIAL CLOSURE: Koide ~ derived from mass formula")
    print("           but residual 0.1-0.5% is not eliminable in current model.")
else:
    print(f"  RESULT:  mass-derived and Koide-derived g_0^tau disagree by {diff_pct:+.3f}%.")
    print("           Koide K=2/3 is NOT simply derived from PDG m_tau/m_e.")
    print("           g_0^tau remains an ADMITTED free parameter.")

print()
print("=" * 78)
print("  ps2 complete.")
print("=" * 78)
