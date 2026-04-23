#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c38_eta_physical_check.py:
Compare perturbation expansion η(δ) against exact ODE for physical leptons.

Known coefficients:
  α₂ = c₁/2 = (1-ln3/4)/2 = 0.362673... (PROVEN)
  α₃ = π²/128 + P_cos = 0.089722... (30 digits)
  α₄ ≈ -0.02460 (5 digits)
  α₅ ≈ 0.02751 (5 digits)

Physical values (from Compton matching + phi-ladder):
  g₀^e = 0.86941 → δ_e = -0.13059
  g₀^μ = 1.40741 → δ_μ = +0.40741

Test: how well does the truncated series reproduce exact η?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  PERTURBATION SERIES vs EXACT ODE for physical leptons")
print("="*78)

# Known coefficients
ln3 = math.log(3)
alpha2 = (1 - ln3/4) / 2
alpha3 = 0.089722223674
alpha4 = -0.02460
alpha5 = 0.02751
# Rough estimates from earlier work:
alpha6 = 0.0  # unknown
alpha7 = 0.011  # rough
alpha9 = 0.007  # rough

print(f"\n  Perturbation coefficients:")
print(f"    α₂ = {alpha2:.12f} (PROVEN)")
print(f"    α₃ = {alpha3:.12f} (30 digits)")
print(f"    α₄ = {alpha4:.6f} (5 digits)")
print(f"    α₅ = {alpha5:.6f} (5 digits)")

def eta_series(delta, order=5):
    """Perturbation series η(δ) truncated at given order."""
    d = delta
    eta = 1.0
    if order >= 2: eta += alpha2 * d
    if order >= 3: eta += alpha3 * d**2
    if order >= 4: eta += alpha4 * d**3
    if order >= 5: eta += alpha5 * d**4
    return eta

# --- Exact ODE ---
def solve_ode(g0, r_max=2000, fit_start=400, fit_end=1900):
    """Solve substrate ODE, extract tail amplitude A_tail and η."""
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]

    r0 = 1e-6
    y0 = [g0, 0.0]

    r_eval = np.linspace(r0, r_max, r_max * 10 + 1)
    sol = solve_ivp(rhs, [r0, r_max], y0, method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)

    if not sol.success:
        return None, None

    r = sol.t
    g = sol.y[0]

    # Tail: g-1 ~ A*sin(r+phi)/r, so r*(g-1) ~ A*sin(r+phi) = B*sin(r) + C*cos(r)
    tail = r * (g - 1)
    mask = (r >= fit_start) & (r <= fit_end)
    r_fit = r[mask]
    t_fit = tail[mask]

    M = np.column_stack([np.sin(r_fit), np.cos(r_fit)])
    coeffs, _, _, _ = np.linalg.lstsq(M, t_fit, rcond=None)
    B, C = coeffs
    A_tail = math.sqrt(B**2 + C**2)

    delta = g0 - 1
    eta = A_tail / abs(delta) if abs(delta) > 1e-10 else 1.0
    return A_tail, eta

# --- Physical lepton parameters ---
phi = (1 + math.sqrt(5)) / 2  # golden ratio
g0_e = 0.86941  # electron
g0_mu = phi * g0_e  # muon (phi-ladder)
g0_tau = 1.72932  # tau (from Koide inversion)

leptons = {
    'electron': g0_e,
    'muon': g0_mu,
    'tau': g0_tau,
}

print(f"\n{'='*78}")
print("  Physical lepton predictions")
print("="*78)

for name, g0 in leptons.items():
    delta = g0 - 1
    print(f"\n  {name}: g₀ = {g0:.5f}, δ = {delta:+.5f}")

    # Exact ODE
    A_tail, eta_exact = solve_ode(g0)
    if eta_exact is None:
        print(f"    ODE failed")
        continue

    print(f"    A_tail (exact) = {A_tail:.8f}")
    print(f"    η (exact)      = {eta_exact:.8f}")

    # Perturbation series at different orders
    for order in [1, 2, 3, 4, 5]:
        eta_pert = eta_series(delta, order)
        err = eta_pert - eta_exact
        rel_err = err / eta_exact * 100
        print(f"    η (order {order})    = {eta_pert:.8f}  err = {err:+.6f} ({rel_err:+.3f}%)")

# --- Scan over delta range ---
print(f"\n{'='*78}")
print("  Series convergence scan")
print("="*78)

delta_values = [0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]

print(f"\n  {'delta':>6s}  {'eta_exact':>12s}  {'O(d)':>12s}  {'O(d2)':>12s}  {'O(d3)':>12s}  {'O(d4)':>12s}  {'err_O4':>10s}")
for delta in delta_values:
    g0 = 1 + delta
    _, eta_ex = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta_ex is None: continue

    e1 = eta_series(delta, 2)
    e2 = eta_series(delta, 3)
    e3 = eta_series(delta, 4)
    e4 = eta_series(delta, 5)
    err = e4 - eta_ex

    print(f"  {delta:6.3f}  {eta_ex:12.8f}  {e1:12.8f}  {e2:12.8f}  {e3:12.8f}  {e4:12.8f}  {err:+10.6f}")

# --- eta_sym(delta) check ---
print(f"\n{'='*78}")
print("  eta_sym(delta) = (eta(+d) + eta(-d))/2  vs  1 + α₃d² + α₅d⁴")
print("="*78)

print(f"\n  {'delta':>6s}  {'eta_sym_ODE':>14s}  {'1+a3d2':>14s}  {'1+a3d2+a5d4':>14s}  {'err_O4':>10s}")
for delta in [0.02, 0.05, 0.10, 0.15, 0.20, 0.30]:
    _, eta_plus = solve_ode(1+delta, r_max=1000, fit_start=200, fit_end=900)
    _, eta_minus = solve_ode(1-delta, r_max=1000, fit_start=200, fit_end=900)
    if eta_plus is None or eta_minus is None: continue

    eta_sym = (eta_plus + eta_minus) / 2
    eta_pert2 = 1 + alpha3 * delta**2
    eta_pert4 = 1 + alpha3 * delta**2 + alpha5 * delta**4

    err = eta_pert4 - eta_sym
    print(f"  {delta:6.3f}  {eta_sym:14.10f}  {eta_pert2:14.10f}  {eta_pert4:14.10f}  {err:+10.7f}")

# --- Mass ratio prediction ---
print(f"\n{'='*78}")
print("  MASS RATIO r₂₁ = m_μ/m_e from perturbation series")
print("="*78)

delta_e = g0_e - 1
delta_mu = g0_mu - 1

eta_e_exact = solve_ode(g0_e)[1]
eta_mu_exact = solve_ode(g0_mu)[1]

print(f"  δ_e = {delta_e:.6f}, δ_μ = {delta_mu:.6f}")
print(f"  η_e (exact) = {eta_e_exact:.8f}")
print(f"  η_μ (exact) = {eta_mu_exact:.8f}")

for order in [2, 3, 4, 5]:
    eta_e_pert = eta_series(delta_e, order)
    eta_mu_pert = eta_series(delta_mu, order)

    r21_pert = (delta_mu/delta_e)**4 * (eta_mu_pert/eta_e_pert)**4
    r21_exact = (delta_mu/delta_e)**4 * (eta_mu_exact/eta_e_exact)**4

    print(f"\n  Order {order}:")
    print(f"    η_e = {eta_e_pert:.8f} (err: {eta_e_pert-eta_e_exact:+.6f})")
    print(f"    η_μ = {eta_mu_pert:.8f} (err: {eta_mu_pert-eta_mu_exact:+.6f})")
    print(f"    r₂₁ = {r21_pert:.4f}")

r21_exact = (delta_mu/delta_e)**4 * (eta_mu_exact/eta_e_exact)**4
print(f"\n  r₂₁ (exact ODE) = {r21_exact:.4f}")
print(f"  r₂₁ (PDG m_μ/m_e) = 206.768")

print(f"\n{'='*78}")
