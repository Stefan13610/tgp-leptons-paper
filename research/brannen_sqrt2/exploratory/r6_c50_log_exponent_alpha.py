#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c50_log_exponent_alpha.py:
Deep investigation of the logarithmic singularity exponent.

From r6_c49: η ≈ C·(-ln ε)^α + D with α ≈ 0.390, C ≈ 0.511, D ≈ 0.862
Question: Is α = 2/5 exactly?

Strategy:
1. Much denser grid near δ_crit (ε down to 0.001)
2. Higher integration precision (R_max=3000)
3. Fix α = 2/5, fit only C, D → compare RMSE
4. Fix α = 1/e, 1/3, 3/8 etc → systematic comparison
5. Richardson-like extrapolation: successive fits on shrinking ε windows
6. Test alternative: η = C·(-ln ε)^(2/5) · (1 + a₁/(-ln ε) + ...) + D
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, minimize_scalar
import math

print("="*78)
print("  LOGARITHMIC EXPONENT α: Is α = 2/5?")
print("="*78)

# --- ODE solver with adjustable precision ---
def solve_ode(g0, r_max=3000, fit_start=500, fit_end=2800):
    """Solve TGP ODE, return η = A_tail/|δ|."""
    def rhs(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    sol = solve_ivp(rhs, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=2.3e-14, atol=1e-20, dense_output=True, max_step=0.1)
    if not sol.success: return None

    # Sample for tail fit
    r_fit = np.linspace(fit_start, fit_end, 20001)
    g_fit = sol.sol(r_fit)[0]
    tail = r_fit * (g_fit - 1)

    # 6-term fit: B·sin + A·cos + C·sin/r + D·cos/r + E·sin/r² + F·cos/r²
    M = np.column_stack([
        np.sin(r_fit), np.cos(r_fit),
        np.sin(r_fit)/r_fit, np.cos(r_fit)/r_fit,
        np.sin(r_fit)/r_fit**2, np.cos(r_fit)/r_fit**2
    ])
    coeffs, _, _, _ = np.linalg.lstsq(M, tail, rcond=None)
    A_tail = math.sqrt(coeffs[0]**2 + coeffs[1]**2)
    delta = g0 - 1
    return A_tail / abs(delta)

# --- Collect data with high precision ---
delta_crit = 1.206187

print(f"\n  δ_crit = {delta_crit:.6f}")
print(f"  Collecting η(δ) data with R_max=3000, 6-term fit...")

epsilons_all = sorted(set(
    list(np.arange(0.80, 0.10, -0.02)) +
    list(np.arange(0.10, 0.01, -0.005)) +
    list(np.arange(0.010, 0.002, -0.001)) +
    list(np.arange(0.0050, 0.0010, -0.0005)) +
    [0.0015, 0.0012, 0.0010]
), reverse=True)

deltas = []
etas = []

for eps in epsilons_all:
    delta = delta_crit - eps
    if delta <= 0: continue
    g0 = 1 + delta
    eta = solve_ode(g0)
    if eta is not None and eta > 0:
        deltas.append(delta)
        etas.append(eta)
        if eps <= 0.01 or eps >= 0.5 or abs(eps - 0.1) < 0.001:
            print(f"    ε = {eps:.5f}, δ = {delta:.5f}, η = {eta:.8f}")

deltas = np.array(deltas)
etas = np.array(etas)
epsilons = delta_crit - deltas
print(f"  Collected {len(deltas)} data points")

# ============================================
# PART 1: Fit α on different ε windows
# ============================================
print(f"\n{'='*78}")
print("  PART 1: α from fits on narrowing ε windows")
print("="*78)

def model_log(eps, C, alpha, D):
    return C * (-np.log(eps))**alpha + D

print(f"\n  {'ε_max':>8s}  {'N pts':>6s}  {'α':>10s}  {'C':>10s}  {'D':>10s}  {'RMSE':>10s}")
for eps_max in [0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01]:
    mask = epsilons < eps_max
    if sum(mask) < 5: continue
    try:
        popt, _ = curve_fit(model_log, epsilons[mask], etas[mask],
                           p0=[0.5, 0.4, 0.9], maxfev=10000)
        pred = model_log(epsilons[mask], *popt)
        rmse = np.sqrt(np.mean((etas[mask] - pred)**2))
        print(f"  {eps_max:8.4f}  {sum(mask):6d}  {popt[1]:10.6f}  {popt[0]:10.6f}  {popt[2]:10.6f}  {rmse:10.2e}")
    except:
        print(f"  {eps_max:8.4f}  {sum(mask):6d}  FIT FAILED")

# ============================================
# PART 2: Test specific α values
# ============================================
print(f"\n{'='*78}")
print("  PART 2: Fixed α candidates — fit only C, D")
print("="*78)

mask_narrow = epsilons < 0.1
eps_n = epsilons[mask_narrow]
eta_n = etas[mask_narrow]

candidates = {
    '2/5 = 0.4000': 2/5,
    '3/8 = 0.3750': 3/8,
    '1/e ≈ 0.3679': 1/math.e,
    '0.3902 (free)': 0.3902,
    '1/3 = 0.3333': 1/3,
    'ln2/2≈0.3466': math.log(2)/2,
    '√(3/20)≈0.387': math.sqrt(3/20),
    '(√5-1)/π≈0.394': (math.sqrt(5)-1)/math.pi,
    '2·ln φ≈0.9624/... skip': None,
    'π/8 ≈ 0.3927': math.pi/8,
}

print(f"\n  {'Candidate':>25s}  {'α':>8s}  {'C':>10s}  {'D':>10s}  {'RMSE':>12s}  {'Δ vs free':>10s}")

# First get free fit RMSE
popt_free, _ = curve_fit(model_log, eps_n, eta_n, p0=[0.5, 0.4, 0.9], maxfev=10000)
pred_free = model_log(eps_n, *popt_free)
rmse_free = np.sqrt(np.mean((eta_n - pred_free)**2))
alpha_free = popt_free[1]

print(f"  {'FREE FIT':>25s}  {alpha_free:8.5f}  {popt_free[0]:10.6f}  {popt_free[2]:10.6f}  {rmse_free:12.2e}  {'---':>10s}")

for name, alpha_val in candidates.items():
    if alpha_val is None: continue
    def model_fixed(eps, C, D):
        return C * (-np.log(eps))**alpha_val + D
    try:
        popt_f, _ = curve_fit(model_fixed, eps_n, eta_n, p0=[0.5, 0.9], maxfev=10000)
        pred_f = model_fixed(eps_n, *popt_f)
        rmse_f = np.sqrt(np.mean((eta_n - pred_f)**2))
        ratio = rmse_f / rmse_free
        print(f"  {name:>25s}  {alpha_val:8.5f}  {popt_f[0]:10.6f}  {popt_f[1]:10.6f}  {rmse_f:12.2e}  {ratio:10.4f}x")
    except:
        print(f"  {name:>25s}  {alpha_val:8.5f}  FIT FAILED")

# ============================================
# PART 3: Optimal α via brute force
# ============================================
print(f"\n{'='*78}")
print("  PART 3: Optimal α via grid search")
print("="*78)

def rmse_for_alpha(alpha_val):
    def model_fixed(eps, C, D):
        return C * (-np.log(eps))**alpha_val + D
    try:
        popt_f, _ = curve_fit(model_fixed, eps_n, eta_n, p0=[0.5, 0.9], maxfev=5000)
        pred_f = model_fixed(eps_n, *popt_f)
        return np.sqrt(np.mean((eta_n - pred_f)**2))
    except:
        return 1.0

# Coarse scan
alphas_scan = np.linspace(0.30, 0.50, 201)
rmses_scan = [rmse_for_alpha(a) for a in alphas_scan]
best_idx = np.argmin(rmses_scan)
print(f"\n  Coarse scan [0.30, 0.50]: best α = {alphas_scan[best_idx]:.4f}, RMSE = {rmses_scan[best_idx]:.2e}")

# Fine scan around best
alpha_lo = alphas_scan[best_idx] - 0.01
alpha_hi = alphas_scan[best_idx] + 0.01
alphas_fine = np.linspace(alpha_lo, alpha_hi, 201)
rmses_fine = [rmse_for_alpha(a) for a in alphas_fine]
best_fine = np.argmin(rmses_fine)
alpha_opt = alphas_fine[best_fine]
print(f"  Fine scan [{alpha_lo:.3f}, {alpha_hi:.3f}]: best α = {alpha_opt:.6f}, RMSE = {rmses_fine[best_fine]:.2e}")

# How close is α=2/5?
rmse_2_5 = rmse_for_alpha(0.4)
print(f"\n  α_optimal = {alpha_opt:.6f}")
print(f"  α = 2/5   = 0.400000")
print(f"  Difference: {alpha_opt - 0.4:+.6f}")
print(f"  RMSE(optimal)  = {rmses_fine[best_fine]:.2e}")
print(f"  RMSE(α=2/5)    = {rmse_2_5:.2e}")
print(f"  RMSE ratio     = {rmse_2_5/rmses_fine[best_fine]:.4f}")

# ============================================
# PART 4: Window stability of α (Richardson-like)
# ============================================
print(f"\n{'='*78}")
print("  PART 4: Does α converge to a limit as ε_max → 0?")
print("="*78)

print(f"\n  {'ε_max':>8s}  {'ε_min':>8s}  {'N':>4s}  {'α':>10s}  {'C':>10s}  {'D':>10s}")
for eps_max in [0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.06, 0.05, 0.04, 0.03, 0.02]:
    mask = (epsilons < eps_max) & (epsilons > 0.001)
    if sum(mask) < 5: continue
    try:
        popt, _ = curve_fit(model_log, epsilons[mask], etas[mask],
                           p0=[0.5, 0.4, 0.9], maxfev=10000)
        print(f"  {eps_max:8.4f}  {min(epsilons[mask]):8.5f}  {sum(mask):4d}  {popt[1]:10.6f}  {popt[0]:10.6f}  {popt[2]:10.6f}")
    except:
        pass

# ============================================
# PART 5: Log-log derivative (numerical dη/d(ln ε))
# ============================================
print(f"\n{'='*78}")
print("  PART 5: Local logarithmic exponent α_local(ε)")
print("="*78)

# If η = C·(-ln ε)^α + D, then d(η-D)/d(-ln ε) = C·α·(-ln ε)^(α-1)
# So d ln(η-D)/d ln(-ln ε) = α (constant if model exact)
# Use numerical derivative

# Best D from free fit
D_best = popt_free[2]
log_L = -np.log(epsilons)  # L = -ln ε

# Sort by L
idx_s = np.argsort(log_L)
L_s = log_L[idx_s]
eta_s = etas[idx_s]

# Numerical d ln(η-D)/d ln L
print(f"\n  Using D = {D_best:.6f} from free fit")
print(f"  {'L=-ln ε':>10s}  {'ε':>10s}  {'η-D':>10s}  {'α_local':>10s}")

for i in range(1, len(L_s)-1):
    if eta_s[i] - D_best <= 0: continue
    # Central difference
    dL = L_s[i+1] - L_s[i-1]
    d_ln_eta = np.log(eta_s[i+1] - D_best) - np.log(eta_s[i-1] - D_best)
    d_ln_L = np.log(L_s[i+1]) - np.log(L_s[i-1])
    if abs(d_ln_L) < 1e-10: continue
    alpha_loc = d_ln_eta / d_ln_L
    eps_i = np.exp(-L_s[i])
    if eps_i < 0.2:
        print(f"  {L_s[i]:10.4f}  {eps_i:10.6f}  {eta_s[i]-D_best:10.6f}  {alpha_loc:10.6f}")

# ============================================
# PART 6: Sub-leading corrections
# ============================================
print(f"\n{'='*78}")
print("  PART 6: Sub-leading: η = C·L^α·(1 + a₁/L + a₂/L²) + D")
print("="*78)

def model_subleading(eps, C, alpha, D, a1):
    L = -np.log(eps)
    return C * L**alpha * (1 + a1/L) + D

try:
    popt_sub, _ = curve_fit(model_subleading, eps_n, eta_n,
                           p0=[0.5, 0.4, 0.9, 0.1], maxfev=10000)
    pred_sub = model_subleading(eps_n, *popt_sub)
    rmse_sub = np.sqrt(np.mean((eta_n - pred_sub)**2))
    C_s, alpha_s, D_s, a1_s = popt_sub
    print(f"\n  Free fit with sub-leading correction:")
    print(f"    C = {C_s:.8f}")
    print(f"    α = {alpha_s:.8f}")
    print(f"    D = {D_s:.8f}")
    print(f"    a₁ = {a1_s:.8f}")
    print(f"    RMSE = {rmse_sub:.2e}")
    print(f"  → η ≈ {C_s:.4f}·L^{alpha_s:.4f}·(1 + {a1_s:.4f}/L) + {D_s:.4f}")
    print(f"  vs leading only: RMSE = {rmse_free:.2e}")
    print(f"  Improvement: {rmse_free/rmse_sub:.1f}×")
except Exception as e:
    print(f"  Fit failed: {e}")

# Now fix α = 2/5 with sub-leading
def model_sub_fixed(eps, C, D, a1):
    L = -np.log(eps)
    return C * L**0.4 * (1 + a1/L) + D

try:
    popt_sf, _ = curve_fit(model_sub_fixed, eps_n, eta_n,
                          p0=[0.5, 0.9, 0.1], maxfev=10000)
    pred_sf = model_sub_fixed(eps_n, *popt_sf)
    rmse_sf = np.sqrt(np.mean((eta_n - pred_sf)**2))
    print(f"\n  With α = 2/5 FIXED + sub-leading:")
    print(f"    C = {popt_sf[0]:.8f}")
    print(f"    D = {popt_sf[1]:.8f}")
    print(f"    a₁ = {popt_sf[2]:.8f}")
    print(f"    RMSE = {rmse_sf:.2e}")
    print(f"  vs free α sub-leading: RMSE = {rmse_sub:.2e}")
except Exception as e:
    print(f"  Fit failed: {e}")

# ============================================
# PART 7: Constants check
# ============================================
print(f"\n{'='*78}")
print("  PART 7: Check if C, D are recognizable constants")
print("="*78)

# From free fit on ε < 0.1
C_val = popt_free[0]
D_val = popt_free[2]

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
e = math.e
ln2 = math.log(2)
ln3 = math.log(3)

const_checks_C = {
    '1/2': 0.5,
    'π/6': pi/6,
    '1/√π': 1/math.sqrt(pi),
    '1/e': 1/e,
    'ln φ': math.log(phi),
    '2/π': 2/pi,
    '√(1/4)': 0.5,
    '(3-e)/2': (3-e)/2,
}

const_checks_D = {
    '1-1/e': 1 - 1/e,
    'e/π': e/pi,
    'π/e-1/4': pi/e - 0.25,
    '√(3/4)': math.sqrt(3/4),
    '1-ln(3/4)/2': 1 - math.log(3/4)/2,
    'π²/12': pi**2/12,
    'ln 3 - ln 2': ln3 - ln2,
    'φ/2': phi/2,
    '5/6': 5/6,
    '6/7': 6/7,
    '7/8': 7/8,
}

print(f"\n  C ≈ {C_val:.8f}")
for name, val in const_checks_C.items():
    diff = abs(C_val - val)
    if diff < 0.05:
        print(f"    {name:>15s} = {val:.8f}, diff = {diff:.6f}")

print(f"\n  D ≈ {D_val:.8f}")
for name, val in const_checks_D.items():
    diff = abs(D_val - val)
    if diff < 0.05:
        print(f"    {name:>15s} = {val:.8f}, diff = {diff:.6f}")

print(f"\n{'='*78}")
print("  SUMMARY")
print("="*78)
print(f"\n  Best model: η(δ) ≈ C·(-ln(δ_c - δ))^α + D")
print(f"  α_optimal  = {alpha_opt:.6f}")
print(f"  α = 2/5    = 0.400000")
print(f"  |Δα|       = {abs(alpha_opt - 0.4):.6f}")
print(f"  RMSE ratio (2/5 vs opt) = {rmse_2_5/rmses_fine[best_fine]:.4f}")
if abs(alpha_opt - 0.4) < 0.01:
    print(f"\n  ★ α = 2/5 is CONSISTENT within fitting uncertainty")
elif abs(alpha_opt - 0.4) < 0.03:
    print(f"\n  ★ α = 2/5 is MARGINAL — could be shifted by sub-leading corrections")
else:
    print(f"\n  ✗ α = 2/5 is DISFAVORED")

print(f"\n{'='*78}")
