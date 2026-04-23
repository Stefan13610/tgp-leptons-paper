#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c47_pade_convergence.py:
Analytic structure of η(δ) via Padé approximants and convergence analysis.

MOTIVATION:
η(δ) = 1 + α₂δ + α₃δ² + α₄δ³ + α₅δ⁴ + ...
is the KEY function connecting ODE dynamics to lepton masses.
√m_i ∝ |δ_i|·η(δ_i), and K=2/3 is a constraint on mass ratios.

If we understand the analytic structure of η(δ) — its singularities,
poles, convergence radius — we may understand WHY K=2/3 holds.

PLAN:
1. Construct Padé approximants [2/2], [3/1], [1/3] from known α's
2. Find poles & zeros → singularity structure
3. Compare Padé vs Taylor vs exact ODE on full δ range
4. Estimate convergence radius from coefficient ratios
5. Check: does g₀_crit = 2.250 (δ_crit = 1.250) correspond to a Padé pole?
6. Study: what happens to η at δ = 1 (g₀ = 2, i.e. twice vacuum)?

Known coefficients:
  α₂ = (1-ln3/4)/2  = 0.362673...  (PROVEN)
  α₃ = π²/128+P_cos = 0.089722...  (30 digits)
  α₄ = -0.02460      (5 digits)
  α₅ = 0.02751        (5 digits)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  PADÉ APPROXIMANTS & CONVERGENCE of η(δ)")
print("="*78)

# --- Known coefficients ---
ln3 = math.log(3)
alpha = [0] * 10  # alpha[k] = coefficient of δ^(k-1)
alpha[0] = 1.0            # η = 1 + ...
alpha[1] = (1 - ln3/4)/2  # α₂ = c₁/2
alpha[2] = 0.089722223674 # α₃
alpha[3] = -0.02460        # α₄
alpha[4] = 0.02751         # α₅

print(f"\n  Perturbation coefficients:")
for k in range(5):
    print(f"    α_{k} = {alpha[k]:+.10f}")

# Coefficient ratios for convergence radius estimate
print(f"\n  Coefficient ratios |α_{'{k+1}'}/α_k|:")
for k in range(1, 4):
    ratio = abs(alpha[k+1] / alpha[k])
    R_est = 1.0 / ratio
    print(f"    |α_{k+1}/α_{k}| = {ratio:.6f}  → R_est = {R_est:.3f}")

# Root test
print(f"\n  Root test |α_k|^(1/k):")
for k in range(1, 5):
    root = abs(alpha[k])**(1.0/k)
    R_est = 1.0 / root
    print(f"    |α_{k}|^(1/{k}) = {root:.6f}  → R_est = {R_est:.3f}")

# --- Padé approximants ---
print(f"\n{'='*78}")
print("  PADÉ APPROXIMANTS")
print("="*78)

# η(δ) = 1 + α₁δ + α₂δ² + α₃δ³ + α₄δ⁴
# Use notation: c_k = alpha[k], so η = c₀ + c₁δ + c₂δ² + c₃δ³ + c₄δ⁴
c = alpha[:5]

# [2/2] Padé: P(δ) = (p₀ + p₁δ + p₂δ²)/(1 + q₁δ + q₂δ²)
# Conditions: P(δ) = η(δ) + O(δ⁵) → 5 equations for p₀,p₁,p₂,q₁,q₂
# Expand: (p₀ + p₁δ + p₂δ²) = (1 + q₁δ + q₂δ²)(c₀ + c₁δ + c₂δ² + c₃δ³ + c₄δ⁴) + O(δ⁵)
# δ⁰: p₀ = c₀ → p₀ = 1
# δ¹: p₁ = c₁ + c₀q₁ → p₁ = c₁ + q₁
# δ²: p₂ = c₂ + c₁q₁ + c₀q₂ → p₂ = c₂ + c₁q₁ + q₂
# δ³: 0 = c₃ + c₂q₁ + c₁q₂
# δ⁴: 0 = c₄ + c₃q₁ + c₂q₂

# From δ³ and δ⁴: solve for q₁, q₂
# c₃ + c₂q₁ + c₁q₂ = 0
# c₄ + c₃q₁ + c₂q₂ = 0
A_mat = np.array([[c[2], c[1]], [c[3], c[2]]])
b_vec = np.array([-c[3], -c[4]])
q = np.linalg.solve(A_mat, b_vec)
q1, q2 = q

p0 = c[0]
p1 = c[1] + q1
p2 = c[2] + c[1]*q1 + q2

print(f"\n  [2/2] Padé: η ≈ (p₀ + p₁δ + p₂δ²) / (1 + q₁δ + q₂δ²)")
print(f"    p₀ = {p0:.10f}")
print(f"    p₁ = {p1:.10f}")
print(f"    p₂ = {p2:.10f}")
print(f"    q₁ = {q1:.10f}")
print(f"    q₂ = {q2:.10f}")

# Poles of [2/2] Padé
disc = q1**2 - 4*q2
print(f"\n    Denominator: 1 + {q1:.6f}δ + {q2:.6f}δ²")
print(f"    Discriminant: {disc:.8f}")
if disc >= 0:
    pole1 = (-q1 + math.sqrt(disc)) / (2*q2)
    pole2 = (-q1 - math.sqrt(disc)) / (2*q2)
    print(f"    Poles (real): δ = {pole1:.6f}, {pole2:.6f}")
    print(f"    → g₀ = {1+pole1:.6f}, {1+pole2:.6f}")
else:
    re_pole = -q1 / (2*q2)
    im_pole = math.sqrt(-disc) / (2*abs(q2))
    print(f"    Poles (complex): δ = {re_pole:.6f} ± {im_pole:.6f}i")
    print(f"    |δ_pole| = {math.sqrt(re_pole**2 + im_pole**2):.6f}")
    print(f"    → Convergence radius ≈ {math.sqrt(re_pole**2 + im_pole**2):.4f}")

# Zeros of [2/2] Padé
disc_num = p1**2 - 4*p0*p2
print(f"\n    Numerator: {p0:.6f} + {p1:.6f}δ + {p2:.6f}δ²")
if disc_num >= 0:
    z1 = (-p1 + math.sqrt(disc_num)) / (2*p2)
    z2 = (-p1 - math.sqrt(disc_num)) / (2*p2)
    print(f"    Zeros (real): δ = {z1:.6f}, {z2:.6f}")
else:
    re_z = -p1 / (2*p2)
    im_z = math.sqrt(-disc_num) / (2*abs(p2))
    print(f"    Zeros (complex): δ = {re_z:.6f} ± {im_z:.6f}i")

# [3/1] Padé: (p₀ + p₁δ + p₂δ² + p₃δ³) / (1 + q₁δ)
# δ⁰: p₀ = c₀ = 1
# δ¹: p₁ = c₁ + q₁
# δ²: p₂ = c₂ + c₁q₁
# δ³: p₃ = c₃ + c₂q₁
# δ⁴: 0 = c₄ + c₃q₁ → q₁ = -c₄/c₃
q1_31 = -c[4]/c[3]
p0_31 = 1.0
p1_31 = c[1] + q1_31
p2_31 = c[2] + c[1]*q1_31
p3_31 = c[3] + c[2]*q1_31

print(f"\n  [3/1] Padé: η ≈ (p₀ + p₁δ + p₂δ² + p₃δ³) / (1 + q₁δ)")
print(f"    q₁ = {q1_31:.10f}")
print(f"    Pole at δ = {-1/q1_31:.6f}  (g₀ = {1-1/q1_31:.6f})")
print(f"    p₃ = {p3_31:.10f}")

# [1/3] Padé: (p₀ + p₁δ) / (1 + q₁δ + q₂δ² + q₃δ³)
# δ⁰: p₀ = c₀ = 1
# δ¹: p₁ = c₁ + q₁
# δ²: 0 = c₂ + c₁q₁ + q₂
# δ³: 0 = c₃ + c₂q₁ + c₁q₂ + q₃
# δ⁴: 0 = c₄ + c₃q₁ + c₂q₂ + c₁q₃
# From δ²: q₂ = -c₂ - c₁q₁
# From δ³: q₃ = -c₃ - c₂q₁ - c₁q₂ = -c₃ - c₂q₁ - c₁(-c₂-c₁q₁) = -c₃ - c₂q₁ + c₁c₂ + c₁²q₁
# From δ⁴: c₄ + c₃q₁ + c₂(-c₂ - c₁q₁) + c₁(-c₃ - c₂q₁ + c₁c₂ + c₁²q₁) = 0
# Simplify: c₄ + c₃q₁ - c₂² - c₂c₁q₁ - c₁c₃ - c₁c₂q₁ + c₁²c₂ + c₁³q₁ = 0
# q₁(c₃ - 2c₁c₂ + c₁³) = -c₄ + c₂² + c₁c₃ - c₁²c₂
denom_13 = c[3] - 2*c[1]*c[2] + c[1]**3
numer_13 = -c[4] + c[2]**2 + c[1]*c[3] - c[1]**2*c[2]
q1_13 = numer_13 / denom_13
q2_13 = -c[2] - c[1]*q1_13
q3_13 = -c[3] - c[2]*q1_13 - c[1]*q2_13
p0_13 = 1.0
p1_13 = c[1] + q1_13

print(f"\n  [1/3] Padé: η ≈ (1 + p₁δ) / (1 + q₁δ + q₂δ² + q₃δ³)")
print(f"    p₁ = {p1_13:.10f}")
print(f"    q₁ = {q1_13:.10f}")
print(f"    q₂ = {q2_13:.10f}")
print(f"    q₃ = {q3_13:.10f}")

# Find poles of [1/3] by solving cubic denominator
coeffs_13 = [q3_13, q2_13, q1_13, 1.0]  # q₃δ³ + q₂δ² + q₁δ + 1 = 0
roots_13 = np.roots(coeffs_13)
print(f"    Poles:")
for root in roots_13:
    if abs(root.imag) < 1e-10:
        print(f"      δ = {root.real:.6f} (real), g₀ = {1+root.real:.6f}")
    else:
        print(f"      δ = {root.real:.6f} ± {abs(root.imag):.6f}i, |δ| = {abs(root):.6f}")

# --- ODE solution for comparison ---
print(f"\n{'='*78}")
print("  EXACT ODE η(δ) vs PADÉ vs TAYLOR")
print("="*78)

def solve_ode(g0, r_max=2000, fit_start=400, fit_end=1900):
    def rhs_ode(r, y):
        g, gp = y
        if abs(g) < 1e-15: return [gp, 0]
        return [gp, -(1/g)*gp**2 - (2/r)*gp - g + 1]
    r0 = 1e-6
    r_eval = np.linspace(r0, r_max, r_max*10+1)
    sol = solve_ivp(rhs_ode, [r0, r_max], [g0, 0.0], method='DOP853',
                    rtol=2.3e-14, atol=1e-20, t_eval=r_eval, max_step=0.1)
    if not sol.success: return None
    r = sol.t
    tail = r * (sol.y[0] - 1)
    mask = (r >= fit_start) & (r <= fit_end)
    M = np.column_stack([np.sin(r[mask]), np.cos(r[mask]),
                          np.sin(r[mask])/r[mask], np.cos(r[mask])/r[mask]])
    coeffs, _, _, _ = np.linalg.lstsq(M, tail[mask], rcond=None)
    A_tail = math.sqrt(coeffs[0]**2 + coeffs[1]**2)
    delta = g0 - 1
    eta = A_tail / abs(delta) if abs(delta) > 1e-10 else 1.0
    return eta

def eta_taylor(d):
    return sum(alpha[k]*d**k for k in range(5))

def eta_pade22(d):
    num = p0 + p1*d + p2*d**2
    den = 1 + q1*d + q2*d**2
    return num/den if abs(den) > 1e-15 else float('inf')

def eta_pade31(d):
    num = p0_31 + p1_31*d + p2_31*d**2 + p3_31*d**3
    den = 1 + q1_31*d
    return num/den if abs(den) > 1e-15 else float('inf')

def eta_pade13(d):
    num = p0_13 + p1_13*d
    den = 1 + q1_13*d + q2_13*d**2 + q3_13*d**3
    return num/den if abs(den) > 1e-15 else float('inf')

# Scan δ range
delta_values = [-0.15, -0.13, -0.10, -0.05, 0.05, 0.10, 0.15, 0.20, 0.30,
                0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]

print(f"\n  {'delta':>6s}  {'g0':>6s}  {'eta_ODE':>10s}  {'Taylor':>10s}  {'[2/2]':>10s}  {'[3/1]':>10s}  {'[1/3]':>10s}  {'T_err%':>8s}  {'P22_err%':>8s}")

for delta in delta_values:
    g0 = 1 + delta
    if g0 <= 0.01 or g0 >= 2.24:
        continue

    eta_ode = solve_ode(g0)
    if eta_ode is None:
        continue

    et = eta_taylor(delta)
    ep22 = eta_pade22(delta)
    ep31 = eta_pade31(delta)
    ep13 = eta_pade13(delta)

    err_t = (et - eta_ode) / eta_ode * 100
    err_p22 = (ep22 - eta_ode) / eta_ode * 100

    print(f"  {delta:+6.3f}  {g0:6.3f}  {eta_ode:10.6f}  {et:10.6f}  {ep22:10.6f}  {ep31:10.6f}  {ep13:10.6f}  {err_t:+8.3f}  {err_p22:+8.3f}")

# --- Physical lepton values ---
print(f"\n{'='*78}")
print("  PHYSICAL LEPTON PREDICTIONS")
print("="*78)

phi = (1 + math.sqrt(5)) / 2
g0_e = 0.86941
g0_mu = phi * g0_e
g0_tau = 1.72932

for name, g0 in [('electron', g0_e), ('muon', g0_mu), ('tau', g0_tau)]:
    delta = g0 - 1
    eta_ode = solve_ode(g0)
    et = eta_taylor(delta)
    ep22 = eta_pade22(delta)
    ep31 = eta_pade31(delta)

    A_ode = abs(delta) * eta_ode
    A_t = abs(delta) * et
    A_p22 = abs(delta) * ep22

    print(f"\n  {name}: g₀ = {g0:.5f}, δ = {delta:+.5f}")
    print(f"    η(ODE)    = {eta_ode:.8f}")
    print(f"    η(Taylor) = {et:.8f}  (err: {(et-eta_ode)/eta_ode*100:+.4f}%)")
    print(f"    η([2/2])  = {ep22:.8f}  (err: {(ep22-eta_ode)/eta_ode*100:+.4f}%)")
    print(f"    η([3/1])  = {ep31:.8f}  (err: {(ep31-eta_ode)/eta_ode*100:+.4f}%)")

# Mass ratios
print(f"\n  Mass ratio r₂₁ = m_μ/m_e predictions:")
de = g0_e - 1
dm = g0_mu - 1
dt = g0_tau - 1

for label, eta_func in [("Taylor", eta_taylor), ("[2/2]", eta_pade22), ("[3/1]", eta_pade31)]:
    A_e = abs(de) * eta_func(de)
    A_m = abs(dm) * eta_func(dm)
    r21 = (A_m/A_e)**4
    print(f"    {label:>8s}: r₂₁ = {r21:.4f}")

eta_e_ode = solve_ode(g0_e)
eta_m_ode = solve_ode(g0_mu)
A_e_ode = abs(de) * eta_e_ode
A_m_ode = abs(dm) * eta_m_ode
r21_ode = (A_m_ode/A_e_ode)**4
print(f"    {'ODE':>8s}: r₂₁ = {r21_ode:.4f}")
print(f"    {'PDG':>8s}: r₂₁ = 206.768")

# --- Singularity analysis ---
print(f"\n{'='*78}")
print("  SINGULARITY ANALYSIS")
print("="*78)

# g₀_crit ≈ 2.250 is where the ODE loses the soliton solution
# δ_crit ≈ 1.250
print(f"\n  Known critical point: g₀_crit ≈ 2.250, δ_crit ≈ 1.250")
print(f"  [2/2] Padé poles:")
if disc >= 0:
    print(f"    δ = {pole1:.6f} and {pole2:.6f}")
    closest = min(abs(pole1), abs(pole2))
    print(f"    Closest pole |δ| = {closest:.6f}")
    print(f"    g₀_crit(Padé) = {1+min(pole1, pole2, key=abs):.4f}")
else:
    r_conv = math.sqrt((-q1/(2*q2))**2 + (-disc)/(4*q2**2))
    print(f"    Complex poles at |δ| = {r_conv:.6f}")

# Check ODE near critical point
print(f"\n  ODE η near critical point:")
for delta in [0.8, 0.9, 1.0, 1.05, 1.10, 1.15, 1.20, 1.24]:
    g0 = 1 + delta
    eta_ode = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta_ode is not None:
        ep22 = eta_pade22(delta)
        print(f"    δ = {delta:+.3f} (g₀={g0:.3f}): η_ODE = {eta_ode:.6f}, η_[2/2] = {ep22:.6f}")
    else:
        print(f"    δ = {delta:+.3f} (g₀={g0:.3f}): ODE FAILED")

# --- η(δ) = η(-δ) symmetry analysis ---
print(f"\n{'='*78}")
print("  SYMMETRY: η_sym(δ) and η_asym(δ)")
print("="*78)
print(f"\n  η(δ) = η_sym(δ) + η_asym(δ)")
print(f"  η_sym  = 1 + α₃δ² + α₅δ⁴ + ...  (even part)")
print(f"  η_asym = α₂δ + α₄δ³ + ...        (odd part)")

print(f"\n  Even coefficients (from symmetric part of ODE):")
print(f"    α₃ = {alpha[2]:+.10f}")
print(f"    α₅ = {alpha[4]:+.10f}")
print(f"    ratio α₅/α₃ = {alpha[4]/alpha[2]:+.6f}")

print(f"\n  Odd coefficients (from ln(3)/4 asymmetry):")
print(f"    α₂ = {alpha[1]:+.10f}")
print(f"    α₄ = {alpha[3]:+.10f}")
print(f"    ratio α₄/α₂ = {alpha[3]/alpha[1]:+.6f}")

# The odd part arises from the 1/g term in the ODE (deficit/excess asymmetry)
# Check: is there a pattern in the odd coefficients?
print(f"\n  α₂/α₄ = {alpha[1]/alpha[3]:+.6f}")
print(f"  Is α₄/α₂ = -α₃? → {alpha[3]/alpha[1]:+.8f} vs {-alpha[2]:+.8f} (diff: {alpha[3]/alpha[1]+alpha[2]:.6f})")
print(f"  Is α₄ = -α₂·α₃? → {-alpha[1]*alpha[2]:+.8f} vs {alpha[3]:+.8f} (diff: {alpha[3]+alpha[1]*alpha[2]:.6f})")

# --- Koide constraint in terms of η ---
print(f"\n{'='*78}")
print("  KOIDE K=2/3 as CONSTRAINT on η series")
print("="*78)

# K = (√m₁ + √m₂ + √m₃)² / (3(m₁+m₂+m₃))
# √m_i ∝ A_i = |δ_i|·η(δ_i)
# With φ-ladder: δ_μ = φδ_e + (φ-1), δ_τ = ???

print(f"\n  With φ-ladder g₀^μ = φ·g₀^e:")
print(f"    δ_e = g₀^e - 1 = -0.13059")
print(f"    δ_μ = φ·g₀^e - 1 = +0.40741")
print(f"    δ_μ = φ·δ_e + (φ-1) = {phi*(-0.13059) + (phi-1):+.5f}")

# For K=2/3 in Brannen form: √m_i = a(1 + √2·cos(θ+2πi/3))
# This means: A₁/A₂ = (1+√2·cos(θ))/(1+√2·cos(θ+2π/3))
# And: A₁/A₃ = (1+√2·cos(θ))/(1+√2·cos(θ+4π/3))

# So K=2/3 determines the RATIO of A values.
# With A_i = |δ_i|·η(δ_i), this constrains η.

# Question: for FIXED g₀^e (i.e., fixed δ_e, δ_μ), what δ_τ gives K=2/3?
# And: can we express this constraint purely in terms of the η series?

A_e = abs(de) * eta_e_ode
A_m = abs(dm) * eta_m_ode

# For K=2/3 (Brannen B=√2):
# Find θ from A_m/A_e ratio
from scipy.optimize import brentq

def brannen_ratio(theta, i, j):
    ci = math.cos(theta + 2*math.pi*i/3)
    cj = math.cos(theta + 2*math.pi*j/3)
    return (1 + math.sqrt(2)*ci) / (1 + math.sqrt(2)*cj)

def theta_eq(theta):
    return brannen_ratio(theta, 1, 0) - A_m/A_e

# Find theta
from scipy.optimize import fsolve
theta_solutions = []
for t0 in np.linspace(0, 2*math.pi, 20):
    try:
        res = fsolve(theta_eq, t0, full_output=True)
        if res[2] == 1:
            theta = res[0][0] % (2*math.pi)
            if all(abs(theta - s) > 0.01 for s in theta_solutions):
                theta_solutions.append(theta)
    except:
        pass

for theta in sorted(theta_solutions):
    sq_e = 1 + math.sqrt(2)*math.cos(theta)
    sq_m = 1 + math.sqrt(2)*math.cos(theta + 2*math.pi/3)
    sq_t = 1 + math.sqrt(2)*math.cos(theta + 4*math.pi/3)

    if sq_e <= 0 or sq_m <= 0 or sq_t <= 0:
        continue

    A_tau_pred = A_e * sq_t / sq_e
    # So |δ_τ|·η(δ_τ) = A_tau_pred

    print(f"\n  θ = {theta:.6f} ({theta*180/math.pi:.2f}°)")
    print(f"    √m ratios: {sq_e:.6f} : {sq_m:.6f} : {sq_t:.6f}")
    print(f"    A_τ(Koide) = {A_tau_pred:.8f}")

    # Find δ_τ from A_τ
    def atail_eq(delta_tau):
        g0t = 1 + delta_tau
        if g0t <= 0.01 or g0t >= 2.24:
            return 1.0
        eta_t = solve_ode(g0t, r_max=1000, fit_start=200, fit_end=900)
        if eta_t is None:
            return 1.0
        return abs(delta_tau) * eta_t - A_tau_pred

    try:
        delta_tau = brentq(atail_eq, 0.5, 1.23, xtol=1e-5)
        eta_tau = solve_ode(1 + delta_tau, r_max=1000, fit_start=200, fit_end=900)
        g0_tau_koide = 1 + delta_tau

        print(f"    δ_τ(Koide) = {delta_tau:+.6f}")
        print(f"    g₀^τ(Koide) = {g0_tau_koide:.6f}")
        print(f"    η_τ(Koide) = {eta_tau:.6f}")

        # Compare Taylor vs Padé prediction for η_τ
        et_tau = eta_taylor(delta_tau)
        ep22_tau = eta_pade22(delta_tau)
        print(f"    η_τ(Taylor) = {et_tau:.6f} (err: {(et_tau-eta_tau)/eta_tau*100:+.3f}%)")
        print(f"    η_τ([2/2])  = {ep22_tau:.6f} (err: {(ep22_tau-eta_tau)/eta_tau*100:+.3f}%)")

        # What δ_τ would the Padé predict?
        def atail_pade22(delta_tau):
            return abs(delta_tau) * eta_pade22(delta_tau) - A_tau_pred

        try:
            delta_tau_p22 = brentq(atail_pade22, 0.5, 1.5, xtol=1e-6)
            print(f"    δ_τ([2/2] pred) = {delta_tau_p22:+.6f} (diff from ODE: {delta_tau_p22-delta_tau:+.6f})")
        except:
            print(f"    δ_τ([2/2] pred): brentq failed")

    except Exception as e:
        print(f"    Failed: {e}")

# --- Check: η(δ) as δ → δ_crit ---
print(f"\n{'='*78}")
print("  η(δ) BEHAVIOR NEAR CRITICAL POINT")
print("="*78)

print(f"\n  {'delta':>6s}  {'g0':>6s}  {'eta_ODE':>10s}  {'A_tail':>10s}  {'d*eta':>10s}")
for delta in np.arange(0.1, 1.26, 0.05):
    g0 = 1 + delta
    eta_ode = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta_ode is not None:
        A_tail = delta * eta_ode
        print(f"  {delta:6.3f}  {g0:6.3f}  {eta_ode:10.6f}  {A_tail:10.6f}  {delta*eta_ode:10.6f}")
    else:
        print(f"  {delta:6.3f}  {g0:6.3f}  FAILED")

# A_tail(δ) = δ·η(δ): does it have a maximum?
print(f"\n  A_tail(δ) = |δ|·η(δ): looking for maximum...")
A_max = 0
delta_max = 0
for delta in np.arange(0.01, 1.25, 0.01):
    g0 = 1 + delta
    eta_ode = solve_ode(g0, r_max=1000, fit_start=200, fit_end=900)
    if eta_ode is not None:
        A = delta * eta_ode
        if A > A_max:
            A_max = A
            delta_max = delta

print(f"  A_tail maximum: A = {A_max:.6f} at δ = {delta_max:.3f} (g₀ = {1+delta_max:.3f})")

print(f"\n{'='*78}")
