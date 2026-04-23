#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c49_log_singularity.py:
Test whether η(δ) has a LOGARITHMIC singularity at δ_crit.

From r6_c48: the local power exponent β_local DECREASES as δ → δ_crit:
  ε = δ_c - δ = 0.25: β = 0.117
  ε = δ_c - δ = 0.01: β = 0.046
  ε = δ_c - δ = 0.005: β = 0.040

This is the SIGNATURE of a logarithmic singularity:
  η ~ C · (-ln ε)^α  gives β_apparent = α/(-ln ε) → 0

versus power law:
  η ~ C · ε^(-β)     gives β_apparent → const = β

PLAN:
1. Test model: η = C · (-ln(δ_c - δ))^α + D
2. Test model: η = C · (δ_c - δ)^(-β₀) · (-ln(δ_c - δ))^α
3. Determine which model fits best
4. If logarithmic: extract α and check against known constants
5. Also investigate: δ_τ ≈ (2/3)·δ_crit more precisely
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, curve_fit
import math

print("="*78)
print("  LOGARITHMIC SINGULARITY TEST for η(δ)")
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
    if not sol.success: return None
    r = sol.t
    tail = r * (sol.y[0] - 1)
    mask = (r >= fit_start) & (r <= fit_end)
    M = np.column_stack([np.sin(r[mask]), np.cos(r[mask]),
                          np.sin(r[mask])/r[mask], np.cos(r[mask])/r[mask]])
    coeffs, _, _, _ = np.linalg.lstsq(M, tail[mask], rcond=None)
    A_tail = math.sqrt(coeffs[0]**2 + coeffs[1]**2)
    return A_tail / abs(g0 - 1)

# --- Collect η data near critical point ---
delta_crit = 1.206187  # from r6_c48

print(f"\n  δ_crit = {delta_crit:.6f}")
print(f"  Collecting η(δ) data...")

deltas = []
etas = []

# Dense grid near critical
for epsilon in np.concatenate([
    np.arange(0.80, 0.01, -0.02),
    np.arange(0.10, 0.01, -0.005),
    np.arange(0.010, 0.002, -0.001),
    [0.003, 0.004, 0.005, 0.006, 0.008]
]):
    delta = delta_crit - epsilon
    if delta <= 0: continue
    g0 = 1 + delta
    eta = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta is not None and eta > 0:
        deltas.append(delta)
        etas.append(eta)

# Sort
idx = np.argsort(deltas)
deltas = np.array([deltas[i] for i in idx])
etas = np.array([etas[i] for i in idx])
epsilons = delta_crit - deltas

print(f"  Collected {len(deltas)} data points")

# --- Model 1: Pure logarithmic η = C · (-ln ε)^α + D ---
print(f"\n{'='*78}")
print("  MODEL 1: η = C · (-ln ε)^α + D")
print("="*78)

def model_log(eps, C, alpha, D):
    return C * (-np.log(eps))**alpha + D

# Use points with ε < 0.3
mask_fit = epsilons < 0.3
eps_fit = epsilons[mask_fit]
eta_fit = etas[mask_fit]

try:
    popt, pcov = curve_fit(model_log, eps_fit, eta_fit, p0=[0.5, 0.5, 1.0],
                           maxfev=10000)
    C1, alpha1, D1 = popt
    eta_pred1 = model_log(eps_fit, *popt)
    residuals1 = eta_fit - eta_pred1
    rmse1 = np.sqrt(np.mean(residuals1**2))
    print(f"\n  Fit (ε < 0.3, {len(eps_fit)} points):")
    print(f"    C = {C1:.8f}")
    print(f"    α = {alpha1:.8f}")
    print(f"    D = {D1:.8f}")
    print(f"    RMSE = {rmse1:.2e}")
    print(f"  → η ≈ {C1:.4f} · (-ln ε)^{alpha1:.4f} + {D1:.4f}")
except Exception as e:
    print(f"  Fit failed: {e}")
    C1, alpha1, D1, rmse1 = 0, 0, 0, 1

# --- Model 2: Pure power law η = C · ε^(-β) ---
print(f"\n{'='*78}")
print("  MODEL 2: η = C · ε^(-β)")
print("="*78)

def model_power(eps, C, beta):
    return C * eps**(-beta)

try:
    popt2, pcov2 = curve_fit(model_power, eps_fit, eta_fit, p0=[1.3, 0.07],
                              maxfev=10000)
    C2, beta2 = popt2
    eta_pred2 = model_power(eps_fit, *popt2)
    residuals2 = eta_fit - eta_pred2
    rmse2 = np.sqrt(np.mean(residuals2**2))
    print(f"\n  Fit (ε < 0.3, {len(eps_fit)} points):")
    print(f"    C = {C2:.8f}")
    print(f"    β = {beta2:.8f}")
    print(f"    RMSE = {rmse2:.2e}")
    print(f"  → η ≈ {C2:.4f} · ε^({-beta2:.6f})")
except Exception as e:
    print(f"  Fit failed: {e}")
    C2, beta2, rmse2 = 0, 0, 1

# --- Model 3: Combined η = C · ε^(-β₀) · (-ln ε)^α ---
print(f"\n{'='*78}")
print("  MODEL 3: η = C · ε^(-β₀) · (-ln ε)^α")
print("="*78)

def model_combined(eps, C, beta0, alpha):
    return C * eps**(-beta0) * (-np.log(eps))**alpha

try:
    popt3, pcov3 = curve_fit(model_combined, eps_fit, eta_fit,
                              p0=[1.3, 0.01, 0.3], maxfev=10000)
    C3, beta03, alpha3 = popt3
    eta_pred3 = model_combined(eps_fit, *popt3)
    residuals3 = eta_fit - eta_pred3
    rmse3 = np.sqrt(np.mean(residuals3**2))
    print(f"\n  Fit (ε < 0.3, {len(eps_fit)} points):")
    print(f"    C = {C3:.8f}")
    print(f"    β₀ = {beta03:.8f}")
    print(f"    α = {alpha3:.8f}")
    print(f"    RMSE = {rmse3:.2e}")
    print(f"  → η ≈ {C3:.4f} · ε^({-beta03:.6f}) · (-ln ε)^{alpha3:.4f}")
except Exception as e:
    print(f"  Fit failed: {e}")
    C3, beta03, alpha3, rmse3 = 0, 0, 0, 1

# --- Model 4: η = C₁ + C₂·(-ln ε) (linear in ln ε) ---
print(f"\n{'='*78}")
print("  MODEL 4: η = C₁ + C₂·(-ln ε)  [linear in log]")
print("="*78)

log_eps = -np.log(eps_fit)
coeffs_lin = np.polyfit(log_eps, eta_fit, 1)
eta_pred4 = np.polyval(coeffs_lin, log_eps)
rmse4 = np.sqrt(np.mean((eta_fit - eta_pred4)**2))
print(f"  η ≈ {coeffs_lin[1]:.6f} + {coeffs_lin[0]:.6f}·(-ln ε)")
print(f"  RMSE = {rmse4:.2e}")

# --- Model 5: η = C₁ + C₂·(-ln ε)^(1/2) ---
print(f"\n{'='*78}")
print("  MODEL 5: η = C₁ + C₂·√(-ln ε)")
print("="*78)

sqrt_log = np.sqrt(-np.log(eps_fit))
coeffs_sqrt = np.polyfit(sqrt_log, eta_fit, 1)
eta_pred5 = np.polyval(coeffs_sqrt, sqrt_log)
rmse5 = np.sqrt(np.mean((eta_fit - eta_pred5)**2))
print(f"  η ≈ {coeffs_sqrt[1]:.6f} + {coeffs_sqrt[0]:.6f}·√(-ln ε)")
print(f"  RMSE = {rmse5:.2e}")

# --- Comparison ---
print(f"\n{'='*78}")
print("  MODEL COMPARISON")
print("="*78)
print(f"\n  {'Model':>30s}  {'RMSE':>12s}  {'Formula':>40s}")
print(f"  {'Pure log (-ln ε)^α + D':>30s}  {rmse1:12.2e}  η = {C1:.4f}·(-ln ε)^{alpha1:.4f} + {D1:.4f}")
print(f"  {'Pure power ε^(-β)':>30s}  {rmse2:12.2e}  η = {C2:.4f}·ε^(-{beta2:.6f})")
print(f"  {'Combined ε^(-β)(-ln ε)^α':>30s}  {rmse3:12.2e}  η = {C3:.4f}·ε^(-{beta03:.6f})·(-ln ε)^{alpha3:.4f}")
print(f"  {'Linear in ln ε':>30s}  {rmse4:12.2e}  η = {coeffs_lin[1]:.4f} + {coeffs_lin[0]:.4f}·(-ln ε)")
print(f"  {'√(-ln ε)':>30s}  {rmse5:12.2e}  η = {coeffs_sqrt[1]:.4f} + {coeffs_sqrt[0]:.4f}·√(-ln ε)")

# --- Show residuals for best model ---
best_rmse = min(rmse1, rmse2, rmse3, rmse4, rmse5)
print(f"\n  Best RMSE: {best_rmse:.2e}")

print(f"\n  Detailed comparison (best models):")
print(f"  {'ε':>8s}  {'η_ODE':>10s}  {'log':>10s}  {'power':>10s}  {'combined':>10s}  {'linear':>10s}")
for i in range(len(eps_fit)):
    if eps_fit[i] > 0.2 or eps_fit[i] < 0.004:
        continue
    print(f"  {eps_fit[i]:8.5f}  {eta_fit[i]:10.6f}  {eta_pred1[i]:10.6f}  {eta_pred2[i]:10.6f}  {eta_pred3[i]:10.6f}  {eta_pred4[i]:10.6f}")

# --- δ_τ/δ_crit analysis ---
print(f"\n{'='*78}")
print("  δ_τ/δ_crit = 2/3 ANALYSIS")
print("="*78)

# More precise: find g₀^e where δ_τ/δ_crit = EXACTLY 2/3
phi = (1 + math.sqrt(5)) / 2

def koide_K(a1, a2, a3):
    s = a1 + a2 + a3
    return s**2 / (3*(a1**2 + a2**2 + a3**2))

from scipy.optimize import brentq

def delta_tau_ratio(g0e):
    """Compute δ_τ(Koide)/δ_crit for given g₀^e."""
    g0m = phi * g0e
    A_e = solve_ode(g0e)
    A_m = solve_ode(g0m)
    if A_e is None or A_m is None: return None

    de = g0e - 1
    dm = g0m - 1

    def K_eq(g0t):
        A_t = solve_ode(g0t, r_max=1000, fit_start=200, fit_end=900)
        if A_t is None: return -1
        return koide_K(de*A_e, dm*A_m, (g0t-1)*A_t) - 2/3

    # Wait - A_e is eta, not A_tail. Fix:
    A_e_tail = abs(de) * A_e
    A_m_tail = abs(dm) * A_m

    def K_eq2(g0t):
        eta_t = solve_ode(g0t, r_max=1000, fit_start=200, fit_end=900)
        if eta_t is None: return -1
        A_t = abs(g0t-1) * eta_t
        return koide_K(A_e_tail, A_m_tail, A_t) - 2/3

    try:
        g0t = brentq(K_eq2, 1.5, 1 + delta_crit - 0.01, xtol=1e-6)
        dt = g0t - 1
        return dt / delta_crit
    except:
        return None

# Scan g₀^e to find where ratio = 2/3
print(f"\n  {'g0_e':>8s}  {'delta_tau/delta_crit':>20s}  {'diff from 2/3':>14s}")
g0e_list = np.arange(0.83, 0.91, 0.005)
ratios = []
g0e_vals = []

for g0e in g0e_list:
    r = delta_tau_ratio(g0e)
    if r is not None:
        diff = r - 2/3
        ratios.append(r)
        g0e_vals.append(g0e)
        marker = " ← CLOSEST" if abs(diff) < 0.002 else ""
        print(f"  {g0e:8.5f}  {r:20.8f}  {diff:+14.8f}{marker}")

# Interpolate to find exact g₀^e where ratio = 2/3
if len(ratios) > 3:
    from scipy.interpolate import interp1d
    f_interp = interp1d(g0e_vals, ratios, kind='cubic')
    try:
        g0e_exact = brentq(lambda x: f_interp(x) - 2/3,
                           min(g0e_vals) + 0.01, max(g0e_vals) - 0.01)
        print(f"\n  g₀^e where δ_τ/δ_crit = 2/3: g₀^e = {g0e_exact:.6f}")

        # Check what this g₀^e gives for mass ratios
        g0m_exact = phi * g0e_exact
        eta_e_ex = solve_ode(g0e_exact)
        eta_m_ex = solve_ode(g0m_exact)
        de_ex = g0e_exact - 1
        dm_ex = g0m_exact - 1
        if eta_e_ex and eta_m_ex:
            r21 = (dm_ex/de_ex)**4 * (eta_m_ex/eta_e_ex)**4
            print(f"  At this g₀^e: r₂₁ = {r21:.4f} (PDG: 206.768)")
            print(f"  g₀^μ = {g0m_exact:.6f}, g₀^τ(pred) = {1+2/3*delta_crit:.6f}")
    except:
        print(f"\n  Could not find exact crossing")

# --- Does η have a nice functional form? ---
print(f"\n{'='*78}")
print("  FUNCTIONAL FORM of η(δ)")
print("="*78)

# Test: η(δ) ≈ (1 - (δ/δ_c)^n)^(-1/m) for some n, m?
# This would give η → ∞ as δ → δ_c with specific exponent structure.

def model_critical(delta, n, m, A):
    x = delta / delta_crit
    return A * (1 - x**n)**(-1/m)

# Use all data
try:
    popt_c, _ = curve_fit(model_critical, deltas, etas, p0=[2.0, 10.0, 1.0],
                           maxfev=10000, bounds=([0.5, 1, 0.5], [10, 100, 2]))
    n_c, m_c, A_c = popt_c
    eta_pred_c = model_critical(deltas, *popt_c)
    rmse_c = np.sqrt(np.mean((etas - eta_pred_c)**2))
    print(f"\n  η ≈ {A_c:.4f} · (1 - (δ/δ_c)^{n_c:.4f})^(-1/{m_c:.4f})")
    print(f"  RMSE = {rmse_c:.2e}")

    # Effective exponent near critical: β = n/m
    print(f"  Effective β = n/m = {n_c/m_c:.6f}")
except Exception as e:
    print(f"  Fit failed: {e}")

# Test: η = cosh(A·δ) / cosh(0) = some symmetric function + asymmetric correction
print(f"\n  Alternative: η(δ) ≈ 1/(1 - (δ/δ_c)²)^β + correction?")

# Simple test: η·(1-(δ/δ_c)²)^β should be approximately constant near δ_c
for beta_test in [1/4, 1/3, 1/2, 1]:
    vals = [etas[i] * (1 - (deltas[i]/delta_crit)**2)**beta_test
            for i in range(len(deltas)) if deltas[i] > 0.5]
    if len(vals) > 3:
        cv = np.std(vals) / np.mean(vals) * 100
        print(f"    β={beta_test:.4f}: η·(1-x²)^β mean={np.mean(vals):.4f}, CV={cv:.2f}%")

print(f"\n{'='*78}")
