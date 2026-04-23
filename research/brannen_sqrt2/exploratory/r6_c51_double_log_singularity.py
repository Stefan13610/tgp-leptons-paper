#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c51_double_log_singularity.py:
Test ultra-weak singularity models for η(δ) near δ_crit.

From r6_c50: the exponent α in η ~ (-ln ε)^α is NOT constant.
It DECREASES as ε → 0:
  ε_max = 0.2: α = 0.383
  ε_max = 0.02: α = 0.188
  α_local drops from 0.357 to 0.330 at smallest ε

This is signature of ultra-weak divergence. Testing:
1. η = a + b·ln(ln(1/ε))  [double logarithm]
2. η = a + b·ln(1/ε) + c·(ln(1/ε))²  [quadratic in log]
3. η = a + b/(ln(ε_0/ε))  [reciprocal log — bounded!]
4. η = a·exp(b/ln(1/ε))   [essential singularity type]

Also: more precise δ_crit via very fine bisection with R_max=5000.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, brentq
import math

print("="*78)
print("  DOUBLE-LOGARITHMIC SINGULARITY TEST")
print("="*78)

# --- ODE solver ---
def solve_ode(g0, r_max=3000, fit_start=500, fit_end=2800):
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    sol = solve_ivp(rhs, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=2.3e-14, atol=1e-20, dense_output=True, max_step=0.1)
    if not sol.success: return None
    r_fit = np.linspace(fit_start, fit_end, 20001)
    g_fit = sol.sol(r_fit)[0]
    tail = r_fit * (g_fit - 1)
    M = np.column_stack([
        np.sin(r_fit), np.cos(r_fit),
        np.sin(r_fit)/r_fit, np.cos(r_fit)/r_fit,
        np.sin(r_fit)/r_fit**2, np.cos(r_fit)/r_fit**2
    ])
    coeffs, _, _, _ = np.linalg.lstsq(M, tail, rcond=None)
    A_tail = math.sqrt(coeffs[0]**2 + coeffs[1]**2)
    return A_tail / abs(g0 - 1)

def get_gmin(g0, r_max=200):
    """Get minimum of g(r) to check if solution blows up."""
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    sol = solve_ivp(rhs, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=1e-12, atol=1e-15, dense_output=True, max_step=0.05)
    if not sol.success: return None
    r_dense = np.linspace(r0, min(r_max, sol.t[-1]), 10000)
    g_dense = sol.sol(r_dense)[0]
    return np.min(g_dense)

# ============================================
# PART 0: Refine δ_crit more precisely
# ============================================
print(f"\n{'='*78}")
print("  PART 0: Refine δ_crit with very fine bisection")
print("="*78)

g0_lo = 2.206186
g0_hi = 2.206188
print(f"\n  Starting bisection in [{g0_lo}, {g0_hi}]...")

for i in range(60):
    g0_mid = (g0_lo + g0_hi) / 2
    gmin = get_gmin(g0_mid, r_max=300)
    if gmin is None or gmin < 1e-6:
        g0_hi = g0_mid
    else:
        g0_lo = g0_mid
    if i % 10 == 0 or i >= 50:
        print(f"  Step {i:3d}: g₀ ∈ [{g0_lo:.15f}, {g0_hi:.15f}], width = {g0_hi-g0_lo:.2e}")

delta_crit = (g0_lo + g0_hi) / 2 - 1
print(f"\n  δ_crit = {delta_crit:.15f}")
print(f"  g₀_crit = {delta_crit + 1:.15f}")

# ============================================
# PART 1: Collect dense data near δ_crit
# ============================================
print(f"\n{'='*78}")
print("  PART 1: Dense data collection")
print("="*78)

eps_list = sorted(set(
    list(np.logspace(np.log10(0.001), np.log10(0.8), 60)) +
    list(np.logspace(np.log10(0.0005), np.log10(0.003), 15)) +
    [0.0003, 0.0004, 0.0006, 0.0008]
), reverse=True)

deltas = []
etas = []
eps_good = []

for eps in eps_list:
    delta = delta_crit - eps
    if delta <= 0: continue
    g0 = 1 + delta
    eta = solve_ode(g0)
    if eta is not None and eta > 0:
        deltas.append(delta)
        etas.append(eta)
        eps_good.append(eps)

deltas = np.array(deltas)
etas = np.array(etas)
eps_arr = np.array(eps_good)
L = -np.log(eps_arr)  # L = ln(1/ε)

print(f"  Collected {len(deltas)} data points")
print(f"  ε range: [{min(eps_arr):.6f}, {max(eps_arr):.4f}]")
print(f"  L range: [{min(L):.4f}, {max(L):.4f}]")
print(f"  η range: [{min(etas):.6f}, {max(etas):.6f}]")

# ============================================
# PART 2: Model comparison
# ============================================
print(f"\n{'='*78}")
print("  PART 2: Ultra-weak singularity models")
print("="*78)

# Use all data
mask = eps_arr < 0.5  # exclude very far from critical

results = {}

# Model A: η = a + b·L^α  (general power of log)
def mA(eps, a, b, alpha):
    return a + b * (-np.log(eps))**alpha
try:
    p, _ = curve_fit(mA, eps_arr[mask], etas[mask], p0=[0.8, 0.5, 0.4], maxfev=20000)
    pred = mA(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['A: a+b·L^α'] = (rmse, f'α={p[2]:.6f}, a={p[0]:.6f}, b={p[1]:.6f}')
    print(f"  A: η = {p[0]:.6f} + {p[1]:.6f}·L^{p[2]:.6f}  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  A: failed ({e})")

# Model B: η = a + b·ln(L) [double log]
def mB(eps, a, b):
    return a + b * np.log(-np.log(eps))
try:
    p, _ = curve_fit(mB, eps_arr[mask], etas[mask], p0=[1.0, 0.3], maxfev=10000)
    pred = mB(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['B: a+b·ln(L)'] = (rmse, f'a={p[0]:.6f}, b={p[1]:.6f}')
    print(f"  B: η = {p[0]:.6f} + {p[1]:.6f}·ln(L)  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  B: failed ({e})")

# Model C: η = a + b·L + c·L² [quadratic in L]
def mC(eps, a, b, c):
    Lv = -np.log(eps)
    return a + b*Lv + c*Lv**2
try:
    p, _ = curve_fit(mC, eps_arr[mask], etas[mask], p0=[1.0, 0.1, 0.001], maxfev=10000)
    pred = mC(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['C: a+b·L+c·L²'] = (rmse, f'a={p[0]:.6f}, b={p[1]:.6f}, c={p[2]:.6f}')
    print(f"  C: η = {p[0]:.6f} + {p[1]:.6f}·L + {p[2]:.6f}·L²  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  C: failed ({e})")

# Model D: η = a + b·√L [sqrt of log]
def mD(eps, a, b):
    return a + b * np.sqrt(-np.log(eps))
try:
    p, _ = curve_fit(mD, eps_arr[mask], etas[mask], p0=[1.0, 0.3], maxfev=10000)
    pred = mD(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['D: a+b·√L'] = (rmse, f'a={p[0]:.6f}, b={p[1]:.6f}')
    print(f"  D: η = {p[0]:.6f} + {p[1]:.6f}·√L  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  D: failed ({e})")

# Model E: η = a + b·L^(1/3) [cube root of log]
def mE(eps, a, b):
    return a + b * (-np.log(eps))**(1/3)
try:
    p, _ = curve_fit(mE, eps_arr[mask], etas[mask], p0=[0.8, 0.5], maxfev=10000)
    pred = mE(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['E: a+b·L^(1/3)'] = (rmse, f'a={p[0]:.6f}, b={p[1]:.6f}')
    print(f"  E: η = {p[0]:.6f} + {p[1]:.6f}·L^(1/3)  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  E: failed ({e})")

# Model F: η = a + b·ln(L) + c·(ln(L))² [quadratic in double-log]
def mF(eps, a, b, c):
    lnL = np.log(-np.log(eps))
    return a + b*lnL + c*lnL**2
try:
    p, _ = curve_fit(mF, eps_arr[mask], etas[mask], p0=[1.0, 0.3, 0.01], maxfev=10000)
    pred = mF(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['F: a+b·ln(L)+c·ln²(L)'] = (rmse, f'a={p[0]:.6f}, b={p[1]:.6f}, c={p[2]:.6f}')
    print(f"  F: η = {p[0]:.6f} + {p[1]:.6f}·ln(L) + {p[2]:.6f}·ln²(L)  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  F: failed ({e})")

# Model G: η = a + b·L^α + c·L^(α-1) [power of log + sub-leading]
def mG(eps, a, b, alpha, c):
    Lv = -np.log(eps)
    return a + b * Lv**alpha + c * Lv**(alpha-1)
try:
    p, _ = curve_fit(mG, eps_arr[mask], etas[mask], p0=[0.8, 0.5, 0.35, 0.1], maxfev=20000)
    pred = mG(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['G: a+b·L^α+c·L^(α-1)'] = (rmse, f'α={p[2]:.6f}, a={p[0]:.6f}, b={p[1]:.6f}, c={p[3]:.6f}')
    print(f"  G: η = {p[0]:.6f} + {p[1]:.6f}·L^{p[2]:.6f} + {p[3]:.6f}·L^{p[2]-1:.6f}  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  G: failed ({e})")

# Model H: η = a·exp(b·L^c)  [stretched exponential in log]
def mH(eps, a, b, c):
    Lv = -np.log(eps)
    return a * np.exp(b * Lv**c)
try:
    p, _ = curve_fit(mH, eps_arr[mask], etas[mask], p0=[1.0, 0.1, 0.3],
                    maxfev=20000, bounds=([0.5, 0, 0], [2, 1, 1]))
    pred = mH(eps_arr[mask], *p)
    rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
    results['H: a·exp(b·L^c)'] = (rmse, f'a={p[0]:.6f}, b={p[1]:.6f}, c={p[2]:.6f}')
    print(f"  H: η = {p[0]:.6f}·exp({p[1]:.6f}·L^{p[2]:.6f})  RMSE={rmse:.2e}")
except Exception as e:
    print(f"  H: failed ({e})")

# --- Ranking ---
print(f"\n{'='*78}")
print("  MODEL RANKING")
print("="*78)
sorted_models = sorted(results.items(), key=lambda x: x[1][0])
for rank, (name, (rmse, info)) in enumerate(sorted_models, 1):
    marker = " ← BEST" if rank == 1 else ""
    print(f"  {rank}. {name:>30s}  RMSE={rmse:.2e}  {info}{marker}")

# ============================================
# PART 3: Focus on best model — zoom analysis
# ============================================
print(f"\n{'='*78}")
print("  PART 3: Stability of best model across ε windows")
print("="*78)

# Refit best model on increasingly narrow windows
best_name = sorted_models[0][0]
print(f"\n  Best model: {best_name}")

# For each window, fit the general L^α model and track α
print(f"\n  {'Window':>20s}  {'N':>4s}  {'α':>10s}  {'RMSE':>10s}")
for eps_hi, eps_lo in [(0.5, 0.001), (0.3, 0.001), (0.2, 0.001), (0.1, 0.001),
                       (0.05, 0.001), (0.03, 0.001), (0.02, 0.001), (0.01, 0.001),
                       (0.1, 0.01), (0.1, 0.005), (0.05, 0.001), (0.02, 0.002)]:
    mask_w = (eps_arr < eps_hi) & (eps_arr > eps_lo)
    if sum(mask_w) < 5: continue
    try:
        p, _ = curve_fit(mA, eps_arr[mask_w], etas[mask_w], p0=[0.8, 0.5, 0.35], maxfev=20000)
        pred = mA(eps_arr[mask_w], *p)
        rmse = np.sqrt(np.mean((etas[mask_w] - pred)**2))
        print(f"  [{eps_lo:.4f}, {eps_hi:.4f}]  {sum(mask_w):4d}  {p[2]:10.6f}  {rmse:10.2e}")
    except:
        pass

# ============================================
# PART 4: Direct derivative analysis
# ============================================
print(f"\n{'='*78}")
print("  PART 4: d(η)/d(L) vs L — what functional form?")
print("="*78)

# Sort by L (ascending)
idx = np.argsort(L)
L_s = L[idx]
eta_s = etas[idx]

# Numerical derivatives
deta_dL = np.gradient(eta_s, L_s)

print(f"\n  {'L':>8s}  {'ε':>10s}  {'η':>10s}  {'dη/dL':>12s}  {'L·dη/dL':>10s}  {'L²·dη/dL':>10s}")
for i in range(0, len(L_s), max(1, len(L_s)//25)):
    if L_s[i] < 0.5: continue
    eps_i = np.exp(-L_s[i])
    print(f"  {L_s[i]:8.4f}  {eps_i:10.6f}  {eta_s[i]:10.6f}  {deta_dL[i]:12.6f}  {L_s[i]*deta_dL[i]:10.6f}  {L_s[i]**2*deta_dL[i]:10.6f}")

# If η ~ a + b·L^α: dη/dL = b·α·L^(α-1)
# So ln(dη/dL) vs ln(L) should be linear with slope (α-1)
mask_deriv = (L_s > 1) & (deta_dL > 0)
if sum(mask_deriv) > 5:
    ln_deta = np.log(deta_dL[mask_deriv])
    ln_L = np.log(L_s[mask_deriv])
    slope, intercept = np.polyfit(ln_L, ln_deta, 1)
    alpha_from_deriv = slope + 1
    print(f"\n  From ln(dη/dL) vs ln(L): slope = {slope:.6f}")
    print(f"  → α = slope + 1 = {alpha_from_deriv:.6f}")
    print(f"  Deviation from integer/fraction:")
    for name, val in [('1/3', 1/3), ('1/4', 1/4), ('2/5', 2/5), ('0', 0)]:
        print(f"    α - {name} = {alpha_from_deriv - val:+.6f}")

# If η ~ a + b·ln(L): dη/dL = b/L
# So L·dη/dL should be constant
product = L_s[mask_deriv] * deta_dL[mask_deriv]
cv_product = np.std(product) / np.mean(product) * 100
print(f"\n  Test η ~ a + b·ln(L):")
print(f"    L·dη/dL mean = {np.mean(product):.6f}, CV = {cv_product:.1f}%")
print(f"    (If exact, CV would be 0%)")

# ============================================
# PART 5: Is there a FINITE limit?
# ============================================
print(f"\n{'='*78}")
print("  PART 5: Does η actually diverge? Check growth rate")
print("="*78)

# η at smallest ε
print(f"\n  ε → 0 approach:")
print(f"  {'ε':>10s}  {'η':>10s}  {'Δη/Δ(lnε)':>12s}")
for i in range(len(eps_arr)-1, -1, -1):
    if eps_arr[i] < 0.01:
        if i < len(eps_arr)-1:
            deta = etas[i] - etas[i+1]
            dlneps = np.log(eps_arr[i]) - np.log(eps_arr[i+1])
            rate = deta / dlneps if abs(dlneps) > 1e-10 else 0
            print(f"  {eps_arr[i]:10.6f}  {etas[i]:10.6f}  {rate:12.6f}")

# Extrapolate: if η = a + b·L^α with α = 0.35, what's η at ε = 1e-100?
alpha_est = 0.35
p_est, _ = curve_fit(mA, eps_arr[mask], etas[mask], p0=[0.8, 0.5, 0.35], maxfev=20000)
for log_eps in [10, 20, 50, 100, 1000]:
    eta_extrap = p_est[0] + p_est[1] * log_eps**p_est[2]
    print(f"  ε = 10^(-{log_eps:4d}):  η_extrap = {eta_extrap:.4f}  (L = {log_eps*np.log(10):.1f})")

print(f"\n  Even at ε = 10^(-1000), η ≈ {p_est[0] + p_est[1] * (1000*np.log(10))**p_est[2]:.2f}")
print(f"  This is very weak divergence — almost bounded!")

print(f"\n{'='*78}")
print("  SUMMARY")
print("="*78)
print(f"\n  δ_crit = {delta_crit:.15f}")
print(f"\n  Key finding: the singularity exponent α in η ~ L^α")
print(f"  is NOT constant — it DECREASES as ε → 0.")
print(f"  This means the singularity is WEAKER than any L^α power.")
print(f"\n  Best candidates (by RMSE):")
for rank, (name, (rmse, info)) in enumerate(sorted_models[:3], 1):
    print(f"    {rank}. {name}: RMSE={rmse:.2e}")

print(f"\n{'='*78}")
