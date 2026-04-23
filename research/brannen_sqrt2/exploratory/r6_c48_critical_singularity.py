#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c48_critical_singularity.py:
Nature of the singularity in η(δ) at the critical soliton δ = δ_crit.

MOTIVATION:
From r6_c47, η(δ) diverges as δ → δ_crit ≈ 1.24 (g₀ → g₀_crit ≈ 2.24).
The nature of this singularity — power law, logarithmic, essential —
determines the analytic structure of η and may connect to Koide K=2/3.

PLAN:
1. Precisely locate δ_crit (g₀_crit) via bisection
2. Determine singularity type: η ~ (δ_c - δ)^(-β)?
3. Study A_tail(δ) = δ·η(δ) near critical point
4. Check: is δ_crit related to known constants?
5. Study the deficit side: η(δ) for δ → -1 (g₀ → 0)
6. Universal exponents? Compare with mean-field theory

PHYSICAL SIGNIFICANCE:
If η ~ (δ_c - δ)^(-β), then the mass formula √m ∝ δ·η(δ) gives
√m ~ δ·(δ_c - δ)^(-β). The Koide condition constrains three values
of δ on this curve. The exponent β may be calculable from the ODE.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import math

print("="*78)
print("  CRITICAL SINGULARITY of η(δ)")
print("="*78)

# --- ODE solver ---
def solve_ode(g0, r_max=2000, fit_start=400, fit_end=1900):
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    r_eval = np.linspace(r0, r_max, r_max*10+1)
    sol = solve_ivp(rhs, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)
    if not sol.success: return None, None
    r = sol.t
    g = sol.y[0]
    tail = r * (g - 1)
    mask = (r >= fit_start) & (r <= fit_end)
    M = np.column_stack([np.sin(r[mask]), np.cos(r[mask]),
                          np.sin(r[mask])/r[mask], np.cos(r[mask])/r[mask]])
    coeffs, _, _, _ = np.linalg.lstsq(M, tail[mask], rcond=None)
    A_tail = math.sqrt(coeffs[0]**2 + coeffs[1]**2)
    delta = g0 - 1
    eta = A_tail / abs(delta) if abs(delta) > 1e-10 else 1.0
    return A_tail, eta

# --- Step 1: Locate δ_crit precisely ---
print(f"\n  STEP 1: Locate δ_crit via bisection")
print(f"  (g₀_crit = 1 + δ_crit is where the ODE soliton ceases to exist)")

# The soliton exists when g(r) > 0 for all r > 0.
# Near g₀_crit, g(r) dips close to 0 in the core.
# We find g₀_crit as the value where min_r g(r) → 0.

def check_soliton(g0, r_max=500):
    """Returns True if soliton exists (g(r) > 0 for all r)."""
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    r_eval = np.linspace(r0, r_max, r_max*100)
    sol = solve_ivp(rhs, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=1e-12, atol=1e-15, t_eval=r_eval, max_step=0.05)
    if not sol.success: return False
    g_min = np.min(sol.y[0])
    return g_min > 0.01

def find_g_min(g0, r_max=500):
    """Returns minimum of g(r) for r > 0."""
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    r_eval = np.linspace(r0, r_max, r_max*100)
    sol = solve_ivp(rhs, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=1e-12, atol=1e-15, t_eval=r_eval, max_step=0.05)
    if not sol.success: return -1
    # Find minimum in the core region (first half)
    core = sol.y[0][:len(sol.y[0])//2]
    return np.min(core)

# Bisection on g₀ between 2.20 (exists) and 2.30 (doesn't exist)
g_lo, g_hi = 2.20, 2.30
for _ in range(40):
    g_mid = (g_lo + g_hi) / 2
    if check_soliton(g_mid):
        g_lo = g_mid
    else:
        g_hi = g_mid

g0_crit = (g_lo + g_hi) / 2
delta_crit = g0_crit - 1

print(f"  g₀_crit = {g0_crit:.12f}")
print(f"  δ_crit = {delta_crit:.12f}")

# Check g_min near critical
print(f"\n  g_min(g₀) near critical point:")
for g0_test in [g0_crit - 0.01, g0_crit - 0.005, g0_crit - 0.002,
                g0_crit - 0.001, g0_crit - 0.0005, g0_crit]:
    gm = find_g_min(g0_test)
    print(f"    g₀ = {g0_test:.8f}: g_min = {gm:.8f}")

# --- Check if δ_crit is a recognizable constant ---
print(f"\n  δ_crit identification:")
phi = (1 + math.sqrt(5)) / 2
candidates = {
    '5/4': 5/4,
    'φ-1/4': phi - 0.25,
    '9/4-1': 9/4 - 1,
    'φ²/2': phi**2/2,
    'e-1': math.e - 1,
    'π/e': math.pi/math.e,
    '√(3/2)': math.sqrt(3/2),
    'ln(3)': math.log(3),
    'π-2': math.pi - 2,
    '4/π': 4/math.pi,
    '√2': math.sqrt(2),
    'φ-1': phi - 1,
    '(1+√5)/4': (1+math.sqrt(5))/4,
    'π/2.5': math.pi/2.5,
}

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - delta_crit)):
    diff = val - delta_crit
    print(f"    {name:>12s} = {val:.10f}  (diff: {diff:+.8f})")
    if abs(diff) < 0.001:
        print(f"                                         ← CLOSE!")

# Also check g₀_crit
print(f"\n  g₀_crit identification:")
candidates_g0 = {
    '9/4': 9/4,
    'e': math.e,
    '√5': math.sqrt(5),
    '7/3': 7/3,
    '2+1/4': 2.25,
    'φ+φ⁻¹': phi + 1/phi,
    'φ²': phi**2,
    '3-φ': 3 - phi,
    '2+ln(2)/4': 2 + math.log(2)/4,
    '2+1/5': 2.2,
    '2+π/12': 2 + math.pi/12,
}
for name, val in sorted(candidates_g0.items(), key=lambda x: abs(x[1] - g0_crit)):
    diff = val - g0_crit
    print(f"    {name:>14s} = {val:.10f}  (diff: {diff:+.8f})")

# --- Step 2: Singularity exponent ---
print(f"\n{'='*78}")
print("  STEP 2: Singularity exponent β")
print(f"  Testing: η(δ) ~ C·(δ_c - δ)^(-β)")
print("="*78)

# Compute η at many points approaching δ_crit
deltas = []
etas = []
A_tails = []

for g0 in np.arange(1.05, g0_crit - 0.001, 0.01):
    delta = g0 - 1
    A, eta = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta is not None:
        deltas.append(delta)
        etas.append(eta)
        A_tails.append(A)

# Finer grid near critical
for dist in [0.20, 0.15, 0.10, 0.08, 0.06, 0.04, 0.03, 0.02, 0.015, 0.01, 0.008, 0.005]:
    g0 = g0_crit - dist
    delta = g0 - 1
    A, eta = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta is not None:
        deltas.append(delta)
        etas.append(eta)
        A_tails.append(A)

# Sort by delta
idx = np.argsort(deltas)
deltas = [deltas[i] for i in idx]
etas = [etas[i] for i in idx]
A_tails = [A_tails[i] for i in idx]

# Log-log analysis
print(f"\n  {'delta':>8s}  {'dc-d':>10s}  {'eta':>10s}  {'ln(eta)':>10s}  {'ln(dc-d)':>10s}  {'beta_local':>10s}")
for i in range(len(deltas)):
    d = deltas[i]
    dcd = delta_crit - d
    if dcd <= 0: continue
    ln_eta = math.log(etas[i])
    ln_dcd = math.log(dcd)

    # Local exponent from consecutive points
    if i > 0 and delta_crit - deltas[i-1] > 0:
        dcd_prev = delta_crit - deltas[i-1]
        beta_local = -(math.log(etas[i]) - math.log(etas[i-1])) / (math.log(dcd) - math.log(dcd_prev))
    else:
        beta_local = float('nan')

    if dcd < 0.25:
        print(f"  {d:8.5f}  {dcd:10.6f}  {etas[i]:10.6f}  {ln_eta:10.6f}  {ln_dcd:10.6f}  {beta_local:10.6f}")

# Global fit: ln(η) = -β·ln(δ_c - δ) + ln(C)
# Use points with δ_c - δ < 0.3
mask = [(delta_crit - d) < 0.3 and (delta_crit - d) > 0.005 for d in deltas]
x_fit = [math.log(delta_crit - deltas[i]) for i in range(len(deltas)) if mask[i]]
y_fit = [math.log(etas[i]) for i in range(len(deltas)) if mask[i]]

if len(x_fit) > 3:
    coeffs = np.polyfit(x_fit, y_fit, 1)
    beta_global = -coeffs[0]
    C_global = math.exp(coeffs[1])
    print(f"\n  Global power-law fit (δ_c - δ ∈ [0.005, 0.3]):")
    print(f"    β = {beta_global:.8f}")
    print(f"    C = {C_global:.8f}")
    print(f"    η ≈ {C_global:.4f} · (δ_c - δ)^({-beta_global:+.6f})")

    # Is β a recognizable number?
    print(f"\n    β identification:")
    beta_cands = {
        '1/2': 0.5,
        '1/3': 1/3,
        '1/4': 0.25,
        '1/5': 0.2,
        '2/3': 2/3,
        '1/6': 1/6,
        '1/π': 1/math.pi,
        '1/e': 1/math.e,
        'ln2': math.log(2),
        '1/φ': 1/phi,
        '2/π': 2/math.pi,
        '1/√5': 1/math.sqrt(5),
    }
    for name, val in sorted(beta_cands.items(), key=lambda x: abs(x[1] - beta_global)):
        diff = val - beta_global
        print(f"      {name:>6s} = {val:.8f}  (diff: {diff:+.8f})")
        if abs(diff) < 0.01:
            print(f"                                         ← CLOSE!")

# --- Step 3: A_tail = δ·η(δ) near critical ---
print(f"\n{'='*78}")
print("  STEP 3: A_tail(δ) = δ·η(δ) near critical point")
print("="*78)

# If η ~ C(δ_c-δ)^(-β), then A = δ·η ~ Cδ·(δ_c-δ)^(-β)
# A diverges at δ_c if β > 0. The rate is:
# dA/dδ = η + δ·η' = η + δ·(-β)·C·(δ_c-δ)^(-β-1)·(-1) = η + β·δ·η/(δ_c-δ)

print(f"\n  {'delta':>8s}  {'A_tail':>10s}  {'eta':>10s}  {'dA/ddelta':>10s}")
for i in range(len(deltas)):
    d = deltas[i]
    if d < 0.8: continue
    A = A_tails[i]
    if i > 0 and i < len(deltas) - 1:
        dA = (A_tails[i+1] - A_tails[i-1]) / (deltas[i+1] - deltas[i-1])
    else:
        dA = float('nan')
    print(f"  {d:8.5f}  {A:10.6f}  {etas[i]:10.6f}  {dA:10.4f}")

# --- Step 4: Deficit side η(δ) for δ → -1 ---
print(f"\n{'='*78}")
print("  STEP 4: Deficit side η(δ) as δ → -1 (g₀ → 0)")
print("="*78)

print(f"\n  {'delta':>8s}  {'g0':>6s}  {'eta':>10s}  {'A_tail':>10s}")
for delta in [-0.05, -0.10, -0.15, -0.20, -0.30, -0.40, -0.50, -0.60, -0.70,
              -0.75, -0.80, -0.85, -0.90, -0.93, -0.95, -0.97, -0.98, -0.99]:
    g0 = 1 + delta
    A, eta = solve_ode(g0, r_max=2000, fit_start=400, fit_end=1900)
    if eta is not None:
        print(f"  {delta:+8.3f}  {g0:6.3f}  {eta:10.6f}  {A:10.6f}")
    else:
        print(f"  {delta:+8.3f}  {g0:6.3f}  FAILED")

# --- Step 5: Can δ_crit constrain Koide? ---
print(f"\n{'='*78}")
print("  STEP 5: Koide constraint and critical point")
print("="*78)

# Key observation: with φ-ladder and K=2/3, the three δ values are:
# δ_e = g₀^e - 1  (negative, deficit)
# δ_μ = φ·g₀^e - 1 (positive, excess)
# δ_τ = ??? (positive, excess)
#
# The Koide relation K=2/3 constrains δ_τ given δ_e and δ_μ.
# Question: does δ_τ/δ_crit have a nice value?

g0_e = 0.86941
g0_mu = phi * g0_e
de = g0_e - 1
dm = g0_mu - 1

# Solve for δ_τ from Koide
A_e, eta_e = solve_ode(g0_e)
A_m, eta_m = solve_ode(g0_mu)

def koide_K(a1, a2, a3):
    s = a1 + a2 + a3
    return s**2 / (3*(a1**2 + a2**2 + a3**2))

def K_from_g0tau(g0t):
    A_t, _ = solve_ode(g0t, r_max=1000, fit_start=200, fit_end=900)
    if A_t is None: return 0
    return koide_K(A_e, A_m, A_t)

try:
    g0_tau_koide = brentq(lambda g0t: K_from_g0tau(g0t) - 2/3, 1.5, g0_crit - 0.01, xtol=1e-5)
    delta_tau_koide = g0_tau_koide - 1
    K_check = K_from_g0tau(g0_tau_koide)

    print(f"\n  g₀^τ(Koide) = {g0_tau_koide:.8f}")
    print(f"  δ_τ(Koide) = {delta_tau_koide:+.8f}")
    print(f"  K = {K_check:.8f}")
    print(f"  δ_crit = {delta_crit:.8f}")

    # Ratios
    ratio_tau_crit = delta_tau_koide / delta_crit
    ratio_g0_crit = g0_tau_koide / g0_crit
    print(f"\n  Ratios:")
    print(f"    δ_τ/δ_crit = {ratio_tau_crit:.8f}")
    print(f"    g₀^τ/g₀_crit = {ratio_g0_crit:.8f}")
    print(f"    δ_τ/δ_μ = {delta_tau_koide/dm:.8f}")
    print(f"    g₀^τ/g₀^μ = {g0_tau_koide/g0_mu:.8f}")

    # How far is τ from critical?
    margin = delta_crit - delta_tau_koide
    margin_frac = margin / delta_crit
    print(f"\n    Distance from critical: δ_crit - δ_τ = {margin:.6f}")
    print(f"    Fractional: (δ_crit - δ_τ)/δ_crit = {margin_frac:.6f}")

    # Is there a constraint? E.g., δ_τ = (2/3)·δ_crit?
    print(f"\n  Testing: δ_τ = f·δ_crit for simple f:")
    for name, f in [('1/2', 0.5), ('2/3', 2/3), ('3/4', 0.75), ('φ-1', phi-1),
                     ('1/φ', 1/phi), ('√(1/2)', math.sqrt(0.5)), ('ln2', math.log(2)),
                     ('π/4-1/4', (math.pi-1)/4)]:
        pred = f * delta_crit
        diff = pred - delta_tau_koide
        print(f"    {name:>10s}: f·δ_c = {pred:.6f}  (diff: {diff:+.6f})")

except Exception as e:
    print(f"  Failed to find g₀^τ(Koide): {e}")

# --- Step 6: Ratio δ_τ/δ_crit as function of g₀^e ---
print(f"\n{'='*78}")
print("  STEP 6: δ_τ/δ_crit vs g₀^e (universality check)")
print("="*78)

print(f"\n  {'g0_e':>8s}  {'g0_tau':>8s}  {'dt/dc':>8s}  {'gt/gc':>8s}  {'K':>8s}")
for g0e in [0.84, 0.85, 0.86, 0.8694, 0.87, 0.88, 0.89, 0.90]:
    g0m = phi * g0e
    A_e_t, _ = solve_ode(g0e)
    A_m_t, _ = solve_ode(g0m)
    if A_e_t is None or A_m_t is None: continue

    def K_eq(g0t):
        A_t, _ = solve_ode(g0t, r_max=1000, fit_start=200, fit_end=900)
        if A_t is None: return -1
        return koide_K(A_e_t, A_m_t, A_t) - 2/3

    try:
        g0_tau_k = brentq(K_eq, 1.5, g0_crit - 0.01, xtol=1e-5)
        dt_k = g0_tau_k - 1
        K_k = K_eq(g0_tau_k) + 2/3
        print(f"  {g0e:8.5f}  {g0_tau_k:8.5f}  {dt_k/delta_crit:8.5f}  {g0_tau_k/g0_crit:8.5f}  {K_k:8.6f}")
    except:
        print(f"  {g0e:8.5f}  FAILED")

print(f"\n{'='*78}")
