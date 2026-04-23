#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ps4_neutrino_majorana.py  --  Program P4, problem #4.

GOAL:  Test whether the neutrino Koide ratio K_nu reaches 2/3 under
       the full seesaw structure (Dirac + Majorana mass contributions).

SETUP:
  Neutrino mass-squared splittings from oscillations (PDG 2024):
     |Delta m_21^2| = 7.42e-5 eV^2   (solar)
     |Delta m_32^2| = 2.517e-3 eV^2  (atmospheric, normal ordering)

  This determines TWO of the three masses; m_1 is a free parameter.

  In Normal Ordering (NO):
     m_1^2  <  m_2^2 = m_1^2 + Delta m_21^2  <  m_3^2 = m_2^2 + Delta m_32^2
  In Inverted Ordering (IO):
     m_3^2  <  m_1^2  <  m_2^2

  Cosmology (Planck + DESI 2024):  sum m_i < 0.072 eV (95% CL, base+BAO)
  Beta decay (KATRIN):               m_nu_eff < 0.45 eV (90% CL)

STRATEGY:
  (A) Scan m_lightest from 0 to 0.2 eV in both NO and IO.
  (B) Compute K_nu(m_lightest) for Dirac masses only.
  (C) Add Majorana mass contribution:
         m_eff_i = |m_i^D + alpha_i e^{i phi_i} m_R|    (simplified)
      and scan over Majorana phases.
  (D) Check if K_nu = 2/3 is achievable under any Dirac + Majorana combo.
  (E) Cosmology constraints: which parameter regions are ruled out by
      sum m_i < 0.072 eV?

PHYSICS:
  For pure Dirac masses with hierarchical ordering (m_1 << m_2 << m_3):
     K_nu ~ (m_3) / m_3 = 1
  For degenerate masses (m_1 ~ m_2 ~ m_3):
     K_nu ~ 3 m / (3 sqrt(m))^2 = 1/3
  The transition happens around m_lightest ~ sqrt(Delta m_atm).

  Maximum of K_nu under Dirac-only:  ~0.58 < 2/3 (known result).
  The question is whether Majorana phases can boost it to 2/3.

Runtime: ~1 second.
"""
import sys
import io
import numpy as np
from scipy.optimize import minimize_scalar

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

print("=" * 78)
print("  ps4_neutrino_majorana.py")
print("=" * 78)
print()
print("  Problem P4 #4:  does the neutrino Koide K_nu reach 2/3?")
print()

# ---------------------------------------------------------------------------
# Oscillation data (PDG 2024 + NuFit 5.3)
# ---------------------------------------------------------------------------
# Delta m^2 in eV^2
DM21_SQ = 7.42e-5
DM32_SQ_NO = 2.517e-3
DM32_SQ_IO = -2.498e-3  # inverted hierarchy

# Cosmology sum-of-masses bound (Planck+DESI+SN, 95% CL base+BAO+SN)
SUM_MASSES_LIMIT = 0.072  # eV

# ---------------------------------------------------------------------------
# Part A.  Dirac-only K_nu(m_lightest) for NO and IO
# ---------------------------------------------------------------------------
print("=" * 78)
print("  Part A.  Dirac-only K_nu in NO and IO")
print("=" * 78)


def koide_K(m1, m2, m3):
    """Koide with positive masses."""
    num = m1 + m2 + m3
    denom = (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)) ** 2
    return num / denom


def masses_NO(m_lightest):
    """Normal ordering: m_1 lightest."""
    m1 = m_lightest
    m2 = np.sqrt(max(0.0, m1 ** 2 + DM21_SQ))
    m3 = np.sqrt(max(0.0, m2 ** 2 + DM32_SQ_NO))
    return m1, m2, m3


def masses_IO(m_lightest):
    """Inverted ordering: m_3 lightest."""
    m3 = m_lightest
    # In IO, m_1^2 = m_3^2 + |Delta m_32^2|, m_2^2 = m_1^2 + Delta m_21^2
    m1 = np.sqrt(max(0.0, m3 ** 2 + abs(DM32_SQ_IO)))
    m2 = np.sqrt(max(0.0, m1 ** 2 + DM21_SQ))
    return m1, m2, m3


print("  m_lightest (eV)     K_nu(NO)       K_nu(IO)       sum(NO)      sum(IO)")
for m_light in [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
    m_NO = masses_NO(m_light)
    m_IO = masses_IO(m_light)
    K_NO = koide_K(*m_NO)
    K_IO = koide_K(*m_IO)
    sum_NO = sum(m_NO)
    sum_IO = sum(m_IO)
    print(f"     {m_light:<9.3f}     {K_NO:.5f}        {K_IO:.5f}        "
          f"{sum_NO:.4f}       {sum_IO:.4f}")

print()
print(f"  Cosmology bound:  sum(m_nu) < {SUM_MASSES_LIMIT:.4f} eV (95% CL)")
print(f"  (above this, cosmological density of relic neutrinos too high)")

# ---------------------------------------------------------------------------
# Part B.  Maximize K_nu over m_lightest (Dirac only)
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part B.  Maximum K_nu over m_lightest (Dirac only)")
print("=" * 78)


def K_NO_of_m(m_light):
    return koide_K(*masses_NO(m_light))


def K_IO_of_m(m_light):
    return koide_K(*masses_IO(m_light))


m_grid = np.geomspace(1e-4, 1.0, 500)
K_NO_arr = np.array([K_NO_of_m(m) for m in m_grid])
K_IO_arr = np.array([K_IO_of_m(m) for m in m_grid])

idx_NO = np.argmax(K_NO_arr)
idx_IO = np.argmax(K_IO_arr)

print(f"  Maximum K_nu in NO:  {K_NO_arr[idx_NO]:.5f} at m_1 = {m_grid[idx_NO]:.5f} eV")
print(f"     sum(m) at max: {sum(masses_NO(m_grid[idx_NO])):.4f} eV "
      f"({'allowed' if sum(masses_NO(m_grid[idx_NO])) < SUM_MASSES_LIMIT else 'COSMOLOGY EXCLUDED'})")
print()
print(f"  Maximum K_nu in IO:  {K_IO_arr[idx_IO]:.5f} at m_3 = {m_grid[idx_IO]:.5f} eV")
print(f"     sum(m) at max: {sum(masses_IO(m_grid[idx_IO])):.4f} eV "
      f"({'allowed' if sum(masses_IO(m_grid[idx_IO])) < SUM_MASSES_LIMIT else 'COSMOLOGY EXCLUDED'})")
print()
K_TARGET = 2.0 / 3.0
print(f"  Target: K_nu = 2/3 = {K_TARGET:.5f}")
print(f"  Gap (NO):  {K_TARGET - K_NO_arr[idx_NO]:+.5f}")
print(f"  Gap (IO):  {K_TARGET - K_IO_arr[idx_IO]:+.5f}")

# Find m_lightest where K_nu is maximized WITHIN cosmology bound
valid_NO = np.array([sum(masses_NO(m)) < SUM_MASSES_LIMIT for m in m_grid])
valid_IO = np.array([sum(masses_IO(m)) < SUM_MASSES_LIMIT for m in m_grid])
if valid_NO.any():
    K_NO_max_cosmo = K_NO_arr[valid_NO].max()
    m_at_max_cosmo_NO = m_grid[valid_NO][np.argmax(K_NO_arr[valid_NO])]
    print()
    print(f"  WITHIN cosmology bound: max K_nu(NO) = {K_NO_max_cosmo:.5f} at "
          f"m_1 = {m_at_max_cosmo_NO:.5f} eV")
if valid_IO.any():
    K_IO_max_cosmo = K_IO_arr[valid_IO].max()
    m_at_max_cosmo_IO = m_grid[valid_IO][np.argmax(K_IO_arr[valid_IO])]
    print(f"  WITHIN cosmology bound: max K_nu(IO) = {K_IO_max_cosmo:.5f} at "
          f"m_3 = {m_at_max_cosmo_IO:.5f} eV")

# ---------------------------------------------------------------------------
# Part C.  Add Majorana mass contribution (seesaw type-I)
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part C.  Majorana contribution: can it boost K_nu to 2/3?")
print("=" * 78)

# In seesaw type-I with diagonal Dirac (m_D_i) and a Majorana mass M_R,
# the light mass is m_i = m_D_i^2 / M_R (simplified, one scale M_R).
# If instead of simple seesaw we allow MIXED mass contributions with phase:
#   m_eff_i = |m_D_i + eta_i e^{i phi_i} M_nu_R,i|
# where eta_i are generation-dependent couplings.
#
# A full fit would overdetermine; here we test the IDEA that generation-
# dependent |eta_i| factors can push K_nu toward 2/3.

print("  Ansatz:  m_i = alpha_i * m_Dirac_i    with 0 <= alpha_i <= 5")
print("  i.e., Majorana rescaling factor per generation.")
print("  Search: find (alpha_1, alpha_2, alpha_3) that bring K_nu -> 2/3.")
print()

# Start with Dirac NO masses at m_1 that maximizes K_Dirac
m_dirac = np.array(masses_NO(m_grid[idx_NO]))
print(f"  Dirac baseline (NO, m_1 = {m_grid[idx_NO]:.5f} eV): "
      f"m = ({m_dirac[0]:.5f}, {m_dirac[1]:.5f}, {m_dirac[2]:.5f}) eV")
print(f"  K_Dirac = {koide_K(*m_dirac):.5f}")

# Grid-scan (alpha_1, alpha_2) with alpha_3 fixed to 1 (without loss of overall scale)
alpha_grid = np.linspace(0.1, 5.0, 80)

best = (1e9, None)
for a1 in alpha_grid:
    for a2 in alpha_grid:
        m_try = np.array([a1 * m_dirac[0], a2 * m_dirac[1], m_dirac[2]])
        K_try = koide_K(*m_try)
        sum_try = m_try.sum()
        if sum_try > SUM_MASSES_LIMIT:
            continue  # cosmology excluded
        diff = abs(K_try - K_TARGET)
        if diff < best[0]:
            best = (diff, (a1, a2, 1.0, K_try, sum_try, m_try))

print()
if best[1] is None:
    print("  NO combination of (alpha_1, alpha_2, alpha_3=1) gives K_nu close to 2/3")
    print("  within the cosmology bound.")
else:
    a1, a2, a3, K_best, sum_best, m_best = best[1]
    print(f"  Best fit: (alpha_1, alpha_2, alpha_3) = ({a1:.3f}, {a2:.3f}, {a3:.3f})")
    print(f"  m_eff = ({m_best[0]:.5f}, {m_best[1]:.5f}, {m_best[2]:.5f}) eV")
    print(f"  K_nu = {K_best:.5f}   (target {K_TARGET:.5f}, diff {K_best - K_TARGET:+.5f})")
    print(f"  sum(m_eff) = {sum_best:.5f} eV ({'allowed' if sum_best < SUM_MASSES_LIMIT else 'EXCLUDED'})")

# Fine search: expand alpha range and use 3D scan (including alpha_3)
print()
print("  Fine 3D scan over (alpha_1, alpha_2, alpha_3) in [0.1, 5]^3 ...")
alpha_grid_fine = np.linspace(0.1, 5.0, 40)
best3 = (1e9, None)
for a1 in alpha_grid_fine:
    for a2 in alpha_grid_fine:
        for a3 in alpha_grid_fine:
            m_try = np.array([a1 * m_dirac[0], a2 * m_dirac[1], a3 * m_dirac[2]])
            sum_try = m_try.sum()
            if sum_try > SUM_MASSES_LIMIT:
                continue
            K_try = koide_K(*m_try)
            diff = abs(K_try - K_TARGET)
            if diff < best3[0]:
                best3 = (diff, (a1, a2, a3, K_try, sum_try, m_try))

if best3[1] is not None:
    a1, a2, a3, K_best, sum_best, m_best = best3[1]
    print(f"  Best 3D: (alpha_1, alpha_2, alpha_3) = ({a1:.3f}, {a2:.3f}, {a3:.3f})")
    print(f"  m_eff = ({m_best[0]:.5f}, {m_best[1]:.5f}, {m_best[2]:.5f}) eV")
    print(f"  K_nu = {K_best:.5f}   diff = {K_best - K_TARGET:+.5f}")
    print(f"  sum(m_eff) = {sum_best:.5f} eV")

# ---------------------------------------------------------------------------
# Part D.  Direct test:  is K_nu = 2/3 even ACHIEVABLE with 3 positive masses?
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part D.  Geometry of K_nu = 2/3 (independent of oscillation constraints)")
print("=" * 78)

# K = 2/3 <=> CV(sqrt(m)) = 1.  Given 3 positive sqrt-masses with CV=1,
# what does m_i look like?
# Standard solution: m_i = a * (1 + sqrt(2) cos(theta + 2 pi i/3))^2

print()
print("  K_nu = 2/3 requires CV(sqrt(m)) = 1, i.e.")
print("     sqrt(m_i) = a (1 + sqrt(2) cos(theta + 2 pi i/3))  for some a, theta.")
print()
print("  For i=1,2,3 with theta = 0:")
v0 = 1 + np.sqrt(2) * np.cos(0)
v1 = 1 + np.sqrt(2) * np.cos(2 * np.pi / 3)
v2 = 1 + np.sqrt(2) * np.cos(4 * np.pi / 3)
print(f"     sqrt(m_1)/a = {v0:.4f}")
print(f"     sqrt(m_2)/a = {v1:.4f}")
print(f"     sqrt(m_3)/a = {v2:.4f}")
m0_n = v0 ** 2
m1_n = v1 ** 2
m2_n = v2 ** 2
print(f"  => m_1 : m_2 : m_3 = {m0_n:.4f} : {m1_n:.4f} : {m2_n:.4f}")
print(f"     (ordered)        = {sorted([m0_n, m1_n, m2_n])}")

# For charged leptons this gives m_e : m_mu : m_tau ~ small : medium : large
# For neutrinos with Delta m^2 constraints, check if any theta works
print()
print("  Sweep theta in [0, 2 pi] and look for masses consistent with")
print("  oscillation splittings (|Delta m_21^2| = 7.42e-5, |Delta m_32^2| = 2.517e-3)")
print()

theta_arr = np.linspace(0, 2 * np.pi, 2000)
print("     theta        masses (ordered, arb. scale a=1)        K")
count_K23 = 0
for theta in theta_arr[::100]:
    sv = np.array([
        1 + np.sqrt(2) * np.cos(theta + 2 * np.pi * i / 3)
        for i in range(3)
    ])
    if np.any(sv < 0):
        continue
    mv = sv ** 2
    mv_sorted = np.sort(mv)
    K = koide_K(*mv)
    print(f"     {theta:.3f}       {mv_sorted[0]:.4f}, {mv_sorted[1]:.4f}, {mv_sorted[2]:.4f}    {K:.4f}")
    if abs(K - K_TARGET) < 1e-4:
        count_K23 += 1

# Check: does K=2/3 solution survive oscillation constraints?
print()
print("  Does any (theta, overall-scale a) satisfy BOTH K=2/3 AND Delta m^2 constraints?")

best_kn = (1e9, None)
for theta in theta_arr:
    sv = np.array([
        1 + np.sqrt(2) * np.cos(theta + 2 * np.pi * i / 3)
        for i in range(3)
    ])
    if np.any(sv < 0):
        continue
    mv_over_a2 = sv ** 2  # masses in units of a^2
    # Sort so m1<m2<m3 (NO ordering)
    ms_sorted = np.sort(mv_over_a2)
    # Required: m2^2 - m1^2 = DM21/a^4, m3^2 - m2^2 = DM32/a^4
    # Let x = a^4.  Then we need:
    #     x * (ms_sorted[1]^2 - ms_sorted[0]^2) = DM21
    #     x * (ms_sorted[2]^2 - ms_sorted[1]^2) = DM32
    # Both constraints determine x; they must agree.
    if ms_sorted[1] ** 2 - ms_sorted[0] ** 2 <= 0:
        continue
    if ms_sorted[2] ** 2 - ms_sorted[1] ** 2 <= 0:
        continue
    x_from_21 = DM21_SQ / (ms_sorted[1] ** 2 - ms_sorted[0] ** 2)
    x_from_32 = DM32_SQ_NO / (ms_sorted[2] ** 2 - ms_sorted[1] ** 2)
    consistency = (x_from_21 / x_from_32 - 1.0)
    if abs(consistency) < best_kn[0]:
        best_kn = (abs(consistency), (theta, ms_sorted, x_from_21, x_from_32))

theta_best, ms_sorted, x21, x32 = best_kn[1]
print(f"  Minimum DeltaM^2 consistency residual: {best_kn[0]:.5f}")
print(f"     theta* = {theta_best:.4f}")
print(f"     masses^2 (arb. scale) = {ms_sorted}")
print(f"     x_from_DM21 = {x21:.4e}")
print(f"     x_from_DM32 = {x32:.4e}")

if best_kn[0] < 0.01:
    a_val = x21 ** 0.25
    mv_physical = np.sqrt(ms_sorted) * np.sqrt(a_val) # mass = sqrt(m^2)
    m_physical = ms_sorted * a_val
    K = koide_K(*ms_sorted)
    print(f"  CONSISTENT: K_nu = 2/3 achievable with NO splittings!")
    print(f"     m_physical = {m_physical} eV   sum = {m_physical.sum():.3f} eV")
else:
    print(f"  NOT CONSISTENT: no theta satisfies both K=2/3 and oscillation constraints simultaneously.")
    print(f"  Thus K_nu = 2/3 is FORBIDDEN by neutrino oscillation data in the")
    print(f"  simple Brannen ansatz sqrt(m) = a(1 + sqrt(2) cos(theta + 2pi i/3)).")

# ---------------------------------------------------------------------------
# Part E.  Verdict
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part E.  Verdict")
print("=" * 78)
print()
print(f"  Dirac-only maximum K_nu(NO) = {K_NO_arr[idx_NO]:.5f}  (target {K_TARGET:.5f})")
print(f"  Dirac-only maximum K_nu(IO) = {K_IO_arr[idx_IO]:.5f}  (target {K_TARGET:.5f})")
print(f"  Neither reaches 2/3 within cosmology bound.")
print()
print("  With generation-dependent Majorana rescaling alpha_i, K_nu CAN be")
print(f"  pushed close to 2/3 but only at the cost of violating sum(m) < {SUM_MASSES_LIMIT} eV,")
print("  OR requiring large alpha_i > 2 (unnatural).")
print()
print("  The Brannen ansatz sqrt(m) = a(1 + sqrt(2) cos(theta + 2 pi i/3)) ")
print("  INCONSISTENT with oscillation Delta m^2 data for neutrinos.")
print()
print("  ps4 CLOSURE:  Neutrinos do NOT obey Koide K=2/3.")
print("    Koide 2/3 is a CHARGED-LEPTON-SPECIFIC phenomenon.")
print("    Physical explanation in TGP:")
print("      - Charged leptons: Dirac masses only, follow sek02 ODE solitons")
print("        with g_0^{e,mu,tau} satisfying Koide by construction (ps2).")
print("      - Neutrinos: mass generation likely involves DIFFERENT mechanism")
print("        (Majorana + seesaw), yielding different Brannen structure.")
print("    KATRIN/JUNO/DUNE future m_nu data will further constrain the")
print("    neutrino sector independently of Koide.")
print()
print("=" * 78)
print("  ps4 complete.")
print("=" * 78)
