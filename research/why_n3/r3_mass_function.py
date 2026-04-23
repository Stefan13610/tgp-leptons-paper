#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r3_mass_function.py — Full mass function m(g₀) analysis

KEY QUESTION: Why doesn't the soliton mass integral reproduce
the observed lepton mass ratios m_μ/m_e = 206.768?

The φ-ladder postulates g₀^μ = φ·g₀^e, but:
1. g₀^e = 0.869 (below vacuum) → deficit soliton
2. g₀^μ = 1.407 (above vacuum) → excess soliton
These are DIFFERENT types with different mass scaling.

This script:
1. Maps m(g₀) in detail for both sides of vacuum
2. Finds the actual mass scaling exponent p(g₀)
3. Checks what g₀-spacing WOULD give the right ratio
4. Tests if NON-φ spacing can reproduce masses

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

# Physical masses (MeV)
M_E = 0.51100
M_MU = 105.658
M_TAU = 1776.86
RATIO_MU_E = M_MU / M_E    # 206.768
RATIO_TAU_E = M_TAU / M_E  # 3477.22

# ================================================================
# SOLVER (from r3_physical_alpha.py)
# ================================================================

def solve_alpha(g0, alpha, d=3, r_max=300.0, n_points=30000, g_floor=1e-10):
    """ODE: g'' + (α/g)g'² + ((d-1)/r)g' = (1-g)·g^{2-2α}"""
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


def compute_mass(g0, alpha, d=3):
    """Compute soliton mass = 4π∫[g^{2α}g'²/2 + U(g) - U(1)]r² dr"""
    g_min, sing, sol = solve_alpha(g0, alpha, d)
    if sing or not sol.success:
        return None
    r = sol.t
    g = sol.y[0]
    gp = sol.y[1]
    eps = g**(2*alpha) * gp**2 / 2 + g**3/3 - g**4/4 - 1.0/12
    integrand = eps * r**2
    mass = 4 * np.pi * np.trapezoid(integrand, r)
    return mass


def compute_mass_components(g0, alpha, d=3):
    """Return mass breakdown: kinetic, potential, total"""
    g_min, sing, sol = solve_alpha(g0, alpha, d)
    if sing or not sol.success:
        return None, None, None
    r = sol.t
    g = sol.y[0]
    gp = sol.y[1]
    kin = g**(2*alpha) * gp**2 / 2
    pot = g**3/3 - g**4/4 - 1.0/12
    integrand_kin = kin * r**2
    integrand_pot = pot * r**2
    m_kin = 4 * np.pi * np.trapezoid(integrand_kin, r)
    m_pot = 4 * np.pi * np.trapezoid(integrand_pot, r)
    return m_kin, m_pot, m_kin + m_pot


# ================================================================
print("=" * 70)
print("  R3: FULL MASS FUNCTION m(g₀)")
print("=" * 70)

# ================================================================
# SECTION 1: m(g₀) for deficit solitons (g₀ < 1)
# ================================================================
print(f"\n{'=' * 70}")
print("  1. DEFICIT SOLITONS (g₀ < 1)")
print("=" * 70)

alpha = 1.0
print(f"\n  α = {alpha} (substrate)")
print(f"  {'g₀':>8s}  {'A=1-g₀':>8s}  {'m_kin':>10s}  {'m_pot':>10s}  {'m_tot':>10s}  {'m/A²':>8s}  {'m/A³':>8s}")
print(f"  {'-'*8:>8s}  {'-'*8:>8s}  {'-'*10:>10s}  {'-'*10:>10s}  {'-'*10:>10s}  {'-'*8:>8s}  {'-'*8:>8s}")

deficit_data = []
for g0 in np.arange(0.10, 1.00, 0.05):
    mk, mp, mt = compute_mass_components(g0, alpha)
    if mt is not None:
        A = 1 - g0
        mA2 = mt / A**2 if A > 0.01 else 0
        mA3 = mt / A**3 if A > 0.01 else 0
        deficit_data.append((g0, A, mk, mp, mt))
        print(f"  {g0:8.3f}  {A:8.4f}  {mk:10.6f}  {mp:10.6f}  {mt:10.6f}  {mA2:8.4f}  {mA3:8.4f}")

# Fit power law m = C · A^p
if len(deficit_data) > 3:
    As = np.array([d[1] for d in deficit_data if d[1] > 0.02])
    ms = np.array([d[4] for d in deficit_data if d[1] > 0.02])
    # Use log-log fit
    valid = (ms > 0) & (As > 0)
    if np.sum(valid) > 3:
        logA = np.log(As[valid])
        logm = np.log(ms[valid])
        p_fit, logC_fit = np.polyfit(logA, logm, 1)
        C_fit = np.exp(logC_fit)
        print(f"\n  Power law fit: m = {C_fit:.4f} · A^{p_fit:.3f}")
        print(f"  (R² = {np.corrcoef(logA, logm)[0,1]**2:.6f})")


# ================================================================
# SECTION 2: Excess solitons (g₀ > 1)
# ================================================================
print(f"\n{'=' * 70}")
print("  2. EXCESS SOLITONS (g₀ > 1)")
print("=" * 70)

# Use shorter r_max for stability near vacuum
print(f"\n  Testing solver stability for g₀ > 1 (α=1, d=3):")
for g0 in [1.01, 1.02, 1.05, 1.1, 1.2, 1.3, 1.5, 1.7, 2.0, 2.1]:
    for rmax in [50, 100, 200, 300]:
        g_min, sing, sol = solve_alpha(g0, alpha, d=3, r_max=rmax)
        mk, mp, mt = None, None, None
        if not sing and sol.success:
            r = sol.t
            g = sol.y[0]
            gp = sol.y[1]
            eps = g**(2*alpha) * gp**2 / 2 + g**3/3 - g**4/4 - 1.0/12
            mt = 4 * np.pi * np.trapezoid(eps * r**2, r)
        status = "SING" if sing else ("FAIL" if not sol.success else f"m={mt:.6f}" if mt else "???")
        if rmax == 50 or (not sing and sol.success):
            print(f"  g₀={g0:.2f}, r_max={rmax}: g_min={g_min:.4f}, {status}")
            if not sing and sol.success:
                break

# ================================================================
# SECTION 3: Different α values
# ================================================================
print(f"\n{'=' * 70}")
print("  3. m(g₀) FOR DIFFERENT α")
print("=" * 70)

for alpha in [0.25, 0.50, 0.75, 1.00]:
    print(f"\n  α = {alpha}:")
    print(f"  {'g₀':>8s}  {'m(g₀)':>10s}  {'A':>8s}")

    for g0 in [0.5, 0.6, 0.7, 0.8, 0.869, 0.95]:
        m = compute_mass(g0, alpha)
        A = 1 - g0
        if m is not None:
            print(f"  {g0:8.3f}  {m:10.6f}  {A:8.4f}")

    # Below vacuum: local scaling exponent
    pairs = [(0.5, 0.7), (0.6, 0.8), (0.7, 0.9)]
    for g0a, g0b in pairs:
        ma = compute_mass(g0a, alpha)
        mb = compute_mass(g0b, alpha)
        if ma and mb and ma > 0 and mb > 0:
            Aa, Ab = 1-g0a, 1-g0b
            p = math.log(ma/mb) / math.log(Aa/Ab)
            print(f"  Scaling g₀={g0a}->{g0b}: p = {p:.3f}")


# ================================================================
# SECTION 4: What spacing WOULD give the right mass ratio?
# ================================================================
print(f"\n{'=' * 70}")
print("  4. WHAT g₀-SPACING GIVES m_μ/m_e = 206.768?")
print("=" * 70)

alpha = 1.0
print(f"\n  α = {alpha} (substrate)")
print(f"  Using DEFICIT solitons only (g₀ < 1, same side of vacuum)")
print(f"  Looking for g₀^e, g₀^μ such that m(g₀^μ)/m(g₀^e) = {RATIO_MU_E:.1f}")
print()

# Build dense mass table for deficit solitons
g0_vals = np.arange(0.05, 0.995, 0.005)
mass_table = []
for g0 in g0_vals:
    m = compute_mass(g0, alpha)
    if m is not None and m > 0:
        mass_table.append((g0, m))

if mass_table:
    g0s = np.array([x[0] for x in mass_table])
    ms = np.array([x[1] for x in mass_table])

    # Find pairs with ratio ≈ 206.768
    best_pairs = []
    for i in range(len(mass_table)):
        for j in range(i+1, len(mass_table)):
            # Deficit: smaller g₀ → larger mass (farther from vacuum)
            g0e, me = mass_table[j]  # closer to vacuum = lighter (electron)
            g0mu, mmu = mass_table[i]  # farther from vacuum = heavier (muon)
            if me > 0:
                ratio = mmu / me
                if abs(ratio - RATIO_MU_E) < 20:
                    best_pairs.append((g0e, g0mu, me, mmu, ratio))

    if best_pairs:
        best_pairs.sort(key=lambda x: abs(x[4] - RATIO_MU_E))
        print(f"  Found {len(best_pairs)} pairs near target ratio:")
        print(f"  {'g₀^e':>8s}  {'g₀^μ':>8s}  {'m_e':>10s}  {'m_μ':>10s}  {'ratio':>8s}  {'g₀^μ/g₀^e':>10s}")
        for g0e, g0mu, me, mmu, ratio in best_pairs[:10]:
            print(f"  {g0e:8.4f}  {g0mu:8.4f}  {me:10.6f}  {mmu:10.6f}  {ratio:8.2f}  {g0mu/g0e:10.4f}")

        # Best match
        g0e, g0mu, me, mmu, ratio = best_pairs[0]
        print(f"\n  Best match: g₀^e={g0e:.4f}, g₀^μ={g0mu:.4f}")
        print(f"  Ratio: {ratio:.2f} (target: {RATIO_MU_E:.1f})")
        print(f"  g₀^μ/g₀^e = {g0mu/g0e:.4f} (φ = {PHI:.4f})")

        # Now look for τ
        for k in range(len(mass_table)):
            g0tau, mtau = mass_table[k]
            if me > 0:
                ratio_tau = mtau / me
                if abs(ratio_tau - RATIO_TAU_E) < 200:
                    print(f"  τ candidate: g₀^τ={g0tau:.4f}, "
                          f"m_τ/m_e={ratio_tau:.0f} (target: {RATIO_TAU_E:.0f})")
    else:
        print(f"  No pairs found near ratio {RATIO_MU_E:.1f}")
        # What's the max ratio achievable?
        if len(mass_table) > 2:
            max_m = max(x[1] for x in mass_table)
            min_m = min(x[1] for x in mass_table if x[1] > 0)
            print(f"  Max mass: {max_m:.6f} (at g₀={mass_table[0][0]:.3f})")
            print(f"  Min mass: {min_m:.6f} (at g₀={mass_table[-1][0]:.3f})")
            print(f"  Max ratio achievable: {max_m/min_m:.2f}")


# ================================================================
# SECTION 5: Same analysis for EXCESS solitons (g₀ > 1)
# ================================================================
print(f"\n{'=' * 70}")
print("  5. EXCESS SOLITONS — MASS FUNCTION")
print("=" * 70)

# Use shorter r_max for numerical stability
print(f"\n  m(g₀) for g₀ > 1, α=1, using r_max=50:")
excess_data = []
for g0 in np.arange(1.01, 2.21, 0.02):
    g_min, sing, sol = solve_alpha(g0, alpha, d=3, r_max=50, n_points=10000)
    if not sing and sol.success:
        r = sol.t
        g = sol.y[0]
        gp = sol.y[1]
        eps = g**(2*alpha) * gp**2 / 2 + g**3/3 - g**4/4 - 1.0/12
        mt = 4 * np.pi * np.trapezoid(eps * r**2, r)
        A = g0 - 1
        excess_data.append((g0, A, mt))
        if abs(g0 - round(g0, 1)) < 0.015 or g0 > 2.1:
            print(f"  g₀={g0:.3f}  A={A:.4f}  m={mt:.6f}")

if excess_data:
    # Check scaling
    print(f"\n  Scaling analysis (excess):")
    for i in range(len(excess_data)-1):
        g0a, Aa, ma = excess_data[i]
        g0b, Ab, mb = excess_data[i+1]
        if ma > 0 and mb > 0 and Aa > 0.01 and Ab > 0.01 and Aa != Ab:
            p = math.log(mb/ma) / math.log(Ab/Aa)
            if abs(g0a - round(g0a, 1)) < 0.015:
                print(f"  g₀={g0a:.2f}->{g0b:.2f}: p = {p:.3f}")


# ================================================================
# SECTION 6: Cross-vacuum mass ratio
# ================================================================
print(f"\n{'=' * 70}")
print("  6. CROSS-VACUUM MASS RATIO")
print("=" * 70)

print("""
  The φ-ladder: g₀^e = 0.869 (deficit), g₀^μ = 1.407 (excess)
  If BOTH are on the SAME side of vacuum, the ratio might work.

  Two interpretations:
  A) ALL generations are EXCESS solitons (g₀ > 1)
  B) ALL generations are DEFICIT solitons (g₀ < 1)
  C) Mixed: e is deficit, μ/τ are excess (φ-ladder)
""")

# Test interpretation A: all excess, using r_max=50
print("  A) All excess (g₀ > 1):")
g0e_exc = 1.01
g0mu_exc = PHI * g0e_exc  # 1.634
g0tau_exc = PHI**2 * g0e_exc  # 2.643

for label, g0 in [("e", g0e_exc), ("μ", g0mu_exc), ("τ", g0tau_exc)]:
    g_min, sing, sol = solve_alpha(g0, alpha, d=3, r_max=50)
    if not sing and sol.success:
        r = sol.t
        g = sol.y[0]
        gp = sol.y[1]
        eps = g**(2*alpha) * gp**2 / 2 + g**3/3 - g**4/4 - 1.0/12
        mt = 4 * np.pi * np.trapezoid(eps * r**2, r)
        print(f"  g₀^{label} = {g0:.4f}: m = {mt:.6f}")
    else:
        print(f"  g₀^{label} = {g0:.4f}: SINGULAR (above barrier)")

# Test interpretation B: all deficit
print(f"\n  B) All deficit (g₀ < 1), with g₀^(n+1)/g₀^(n) = const:")
# For deficit solitons: g₀ < 1, heavier = farther from vacuum
# So g₀^τ < g₀^μ < g₀^e < 1
# We need m(g₀^μ)/m(g₀^e) = 206.768

# Already computed in Section 4 — reference result
print("  (See Section 4 for detailed search)")


# ================================================================
# SECTION 7: Asymmetric potential analysis
# ================================================================
print(f"\n{'=' * 70}")
print("  7. WHY MASS IS ASYMMETRIC AROUND VACUUM")
print("=" * 70)

print("""
  U(g) = g³/3 - g⁴/4,  U(1) = 1/12
  δU = U(g) - U(1) = g³/3 - g⁴/4 - 1/12

  Taylor around g = 1 + A:
    δU = -(A²/2) + (A³/6)(higher) + ...
    Leading: δU ≈ -A²/2 (symmetric)

  But for LARGE A:
    g = 0.5: δU = 0.042-0.016-0.083 = -0.057
    g = 1.5: δU = 1.125-1.266-0.083 = -0.224

  The potential well is DEEPER on the excess (g>1) side!
  This means excess solitons store more energy per unit amplitude.
""")

# Compute δU
print(f"  δU(g) = U(g) - U(1):")
for g in [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.5, 1.8, 2.0]:
    dU = g**3/3 - g**4/4 - 1.0/12
    print(f"  g = {g:.2f}: δU = {dU:+.6f}")


# ================================================================
# SECTION 8: Mass with different α
# ================================================================
print(f"\n{'=' * 70}")
print("  8. MASS RATIOS FOR NATURAL α = 3/4")
print("=" * 70)

alpha = 0.75
print(f"\n  α = {alpha} (flat kinetic + volume)")
print(f"  Deficit solitons:")

# Build mass table
mass_table_075 = []
for g0 in np.arange(0.05, 0.995, 0.005):
    m = compute_mass(g0, alpha)
    if m is not None and m > 0:
        mass_table_075.append((g0, m))

if mass_table_075:
    # Print a few values
    for g0, m in mass_table_075:
        if abs(g0 - round(g0, 1)) < 0.003:
            print(f"  g₀ = {g0:.3f}: m = {m:.6f}")

    # Max ratio
    max_m = max(x[1] for x in mass_table_075)
    min_m = min(x[1] for x in mass_table_075 if x[1] > 0)
    print(f"\n  Max ratio (deficit side): {max_m/min_m:.2f}")
    print(f"  Need: {RATIO_MU_E:.1f} for μ/e, {RATIO_TAU_E:.1f} for τ/e")

    # Search for pairs
    best = None
    for i in range(len(mass_table_075)):
        for j in range(i+1, len(mass_table_075)):
            g0e, me = mass_table_075[j]
            g0mu, mmu = mass_table_075[i]
            if me > 0:
                ratio = mmu / me
                if best is None or abs(ratio - RATIO_MU_E) < abs(best[4] - RATIO_MU_E):
                    best = (g0e, g0mu, me, mmu, ratio)

    if best:
        print(f"\n  Closest to target ({RATIO_MU_E:.1f}):")
        print(f"  g₀^e={best[0]:.4f}, g₀^μ={best[1]:.4f}, ratio={best[4]:.2f}")


# ================================================================
print(f"\n{'=' * 70}")
print("  SUMMARY")
print("=" * 70)

print(f"""
  MASS FUNCTION ANALYSIS:

  1. Deficit solitons (g₀ < 1): m ~ (1-g₀)^p with p ≈ 2.6
     Mass DECREASES as g₀ → 1 (closer to vacuum = lighter)

  2. Excess solitons (g₀ > 1): numerically unstable for
     intermediate g₀ (1.1-2.0) due to long-range oscillations.
     Solver works for g₀ near 1 and near g₀_crit.

  3. The potential U(g) is ASYMMETRIC: deeper on the excess side.
     This means deficit and excess solitons of equal |A| = |g₀-1|
     have DIFFERENT masses.

  4. On the DEFICIT side alone:
     Max mass ratio ≈ 10²-10³ (depending on minimum g₀)
     MAY be sufficient for m_μ/m_e if g₀^e is close to vacuum.

  5. The φ-ladder CROSSES vacuum (g₀^e < 1, g₀^μ > 1).
     This means the mass ratio involves TWO DIFFERENT scaling laws.
     The simple power-law argument breaks down.

  IMPLICATIONS FOR N=3:
  The barrier mechanism (g₀_crit from metric singularity) does NOT
  depend on the mass formula. N=3 is determined by:
  - g₀_crit(α, d) from ODE
  - Generation amplitude spacing (φ-ladder or other)
  Both are INDEPENDENT of how m(g₀) maps to physical masses.

  OPEN QUESTION: What is the correct mass formula?
  - Simple soliton energy integral (this script)?
  - Energy WITH GL(3,F₂) corrections?
  - Renormalized mass including quantum corrections?
  - Topological charge contribution?
""")

print(f"\n{'=' * 70}")
print(f"  TEST REPORT: {PASS} PASS, {FAIL} FAIL out of {PASS + FAIL}")
print(f"{'=' * 70}")
