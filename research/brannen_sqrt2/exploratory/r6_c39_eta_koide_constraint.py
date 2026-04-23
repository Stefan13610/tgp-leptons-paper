#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c39_eta_koide_constraint.py:
What constraint does K=2/3 impose on η(δ)?

Given the phi-ladder g₀^μ = φ·g₀^e (δ_μ = φ*δ_e + φ-1),
K = 2/3 determines g₀^τ and thus δ_τ.

Question: is there a pattern in η_e, η_μ, η_τ when K=2/3?

The Koide formula in terms of masses: K = (m_e+m_μ+m_τ)²/(3(m_e²+m_μ²+m_τ²))... no wait.
Actually K = Σmᵢ/(Σ√mᵢ)² · 1/3... no.

Actually: K = (√m₁ + √m₂ + √m₃)² / (3(m₁+m₂+m₃))

In the TGP framework, √mᵢ ∝ A_tail(g₀ⁱ) = |δᵢ|·η(δᵢ).

With Brannen parameterization: √mᵢ = a(1 + B·cos(θ + 2πi/3)), and B=√2 ↔ K=2/3.

Key insight: the perturbation expansion gives η(δ) as a universal function.
So the mass ratios are determined by: √mᵢ/√mⱼ = |δᵢ|η(δᵢ) / |δⱼ|η(δⱼ).

Can we find a relationship that η must satisfy for K=2/3?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import math

print("="*78)
print("  KOIDE CONSTRAINT on η(δ)")
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
    M = np.column_stack([np.sin(r[mask]), np.cos(r[mask])])
    coeffs, _, _, _ = np.linalg.lstsq(M, tail[mask], rcond=None)
    return math.sqrt(coeffs[0]**2 + coeffs[1]**2)

# --- Perturbation coefficients ---
ln3 = math.log(3)
alpha2 = (1 - ln3/4) / 2
alpha3 = 0.089722223674
alpha4 = -0.02460
alpha5 = 0.02751

def eta_pert(delta):
    d = delta
    return 1 + alpha2*d + alpha3*d**2 + alpha4*d**3 + alpha5*d**4

# --- Koide formula ---
def koide_K(m1, m2, m3):
    """Koide K-factor."""
    s1 = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return s1**2 / (3*(m1+m2+m3))

# --- Scan: For fixed g₀^e, vary g₀^τ and find K=2/3 ---
phi = (1 + math.sqrt(5)) / 2

print(f"\n  Scanning g₀^e values, finding g₀^τ for K=2/3:")
print(f"  (g₀^μ = φ·g₀^e, A_i = |δ_i|·η(δ_i))")

# Use perturbation series for fast evaluation
def koide_from_g0e_g0tau(g0e, g0tau):
    """Compute K from g0e, g0mu=phi*g0e, g0tau using perturbation eta."""
    g0mu = phi * g0e
    de = g0e - 1
    dm = g0mu - 1
    dt = g0tau - 1

    # √m ∝ A_tail = |δ|·η(δ)
    a_e = abs(de) * eta_pert(de)
    a_m = abs(dm) * eta_pert(dm)
    a_t = abs(dt) * eta_pert(dt)

    # Masses ∝ A²
    m_e = a_e**2
    m_m = a_m**2
    m_t = a_t**2

    return koide_K(m_e, m_m, m_t)

# For each g0e, find g0tau that gives K=2/3
g0e_values = [0.84, 0.85, 0.86, 0.8694, 0.87, 0.88, 0.89, 0.90]

print(f"\n  {'g0_e':>8s}  {'g0_mu':>8s}  {'g0_tau':>8s}  {'delta_tau':>10s}  {'eta_tau':>10s}  {'K':>8s}")
for g0e in g0e_values:
    try:
        # Find g0tau such that K(g0e, phi*g0e, g0tau) = 2/3
        def K_minus_target(g0t):
            return koide_from_g0e_g0tau(g0e, g0t) - 2/3

        # K decreases with g0tau (roughly), search in [1.3, 2.2]
        g0tau = brentq(K_minus_target, 1.3, 2.2)
        K = koide_from_g0e_g0tau(g0e, g0tau)
        dt = g0tau - 1
        et = eta_pert(dt)
        print(f"  {g0e:8.5f}  {phi*g0e:8.5f}  {g0tau:8.5f}  {dt:+10.5f}  {et:10.6f}  {K:8.6f}")
    except Exception as e:
        print(f"  {g0e:8.5f}  FAILED: {e}")

# --- Key question: what is the constraint on η? ---
print(f"\n{'='*78}")
print("  CONSTRAINT ANALYSIS: what condition on η gives K=2/3?")
print("="*78)

# With √m_i ∝ |δ_i|·η(δ_i) and B=√2 (↔ K=2/3):
# √m_i = a(1 + √2·cos(θ+2πi/3))
# So: |δ_i|·η(δ_i) / |δ_j|·η(δ_j) = (1+√2 cos(θ+2πi/3)) / (1+√2 cos(θ+2πj/3))
# This means: for fixed θ, the ratios η(δ_i)/η(δ_j) are determined.

# Find θ from physical values
g0e = 0.86941
g0mu = phi * g0e  # = 1.40673

# Exact ODE
A_e = solve_ode(g0e)
A_mu = solve_ode(g0mu)
de = g0e - 1
dm = g0mu - 1

print(f"\n  Physical values:")
print(f"    g₀^e = {g0e:.5f}, δ_e = {de:+.6f}, A_e = {A_e:.8f}")
print(f"    g₀^μ = {g0mu:.5f}, δ_μ = {dm:+.6f}, A_μ = {A_mu:.8f}")

# From Brannen: √m_i = a(1 + √2·cos(θ+2πi/3))
# Ratio: A_μ/A_e = (1+√2·cos(θ+2π/3)) / (1+√2·cos(θ))
# If we assign e→i=0, μ→i=1:
# A_e = a(1+√2·cos(θ))
# A_μ = a(1+√2·cos(θ+2π/3))
ratio = A_mu / A_e
print(f"    A_μ/A_e = {ratio:.8f}")

# Find θ from this ratio
# (1+√2·cos(θ+2π/3)) / (1+√2·cos(θ)) = ratio
# Using cos(θ+2π/3) = cos(θ)cos(2π/3) - sin(θ)sin(2π/3) = -cos(θ)/2 - √3 sin(θ)/2
from scipy.optimize import fsolve

def theta_eq(theta):
    c0 = math.cos(theta)
    c1 = math.cos(theta + 2*math.pi/3)
    return (1 + math.sqrt(2)*c1) / (1 + math.sqrt(2)*c0) - ratio

# Search for theta
theta_init = [0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
solutions = set()
for t0 in theta_init:
    try:
        t_sol = fsolve(theta_eq, t0, full_output=True)
        if t_sol[2] == 1:  # converged
            theta = t_sol[0][0] % (2*math.pi)
            if all(abs(theta - s) > 0.01 for s in solutions):
                solutions.add(round(theta, 10))
    except:
        pass

print(f"    θ solutions: {sorted(solutions)}")

for theta in sorted(solutions):
    c0 = math.cos(theta)
    c1 = math.cos(theta + 2*math.pi/3)
    c2 = math.cos(theta + 4*math.pi/3)

    # Predicted masses (relative)
    sq_e = 1 + math.sqrt(2)*c0
    sq_m = 1 + math.sqrt(2)*c1
    sq_t = 1 + math.sqrt(2)*c2

    if sq_e < 0 or sq_m < 0 or sq_t < 0:
        continue

    # Predicted A_tau
    A_tau_pred = A_e * sq_t / sq_e
    # So |δ_τ|·η(δ_τ) = A_tau_pred
    # Need to find g₀^τ such that A_tail(g₀^τ) = A_tau_pred

    print(f"\n    θ = {theta:.6f} ({theta*180/math.pi:.2f}°)")
    print(f"    1+√2·cos(θ) = {sq_e:.6f} (∝ √m_e)")
    print(f"    1+√2·cos(θ+2π/3) = {sq_m:.6f} (∝ √m_μ)")
    print(f"    1+√2·cos(θ+4π/3) = {sq_t:.6f} (∝ √m_τ)")
    print(f"    Predicted A_τ/A_e = {sq_t/sq_e:.6f}")
    print(f"    Predicted A_τ = {A_tau_pred:.6f}")

    # Find g0tau from A_tau_pred using bisection
    def atail_minus_target(g0t):
        at = solve_ode(g0t)
        return at - A_tau_pred if at is not None else 1.0

    try:
        g0tau = brentq(atail_minus_target, 1.4, 2.1, xtol=1e-6)
        A_tau_check = solve_ode(g0tau)
        delta_tau = g0tau - 1
        eta_tau = A_tau_check / abs(delta_tau)

        print(f"    g₀^τ (Koide) = {g0tau:.6f}")
        print(f"    δ_τ = {delta_tau:+.6f}")
        print(f"    η_τ = {eta_tau:.6f}")
        print(f"    η_τ (pert) = {eta_pert(delta_tau):.6f}")

        # Check K
        m_e = A_e**2
        m_m = A_mu**2
        m_t = A_tau_check**2
        K = koide_K(m_e, m_m, m_t)
        print(f"    K = {K:.8f} (target: 0.666667)")

        # Eta ratios
        eta_e = A_e / abs(de)
        eta_mu = A_mu / abs(dm)
        print(f"\n    η_e = {eta_e:.8f}")
        print(f"    η_μ = {eta_mu:.8f}")
        print(f"    η_τ = {eta_tau:.8f}")
        print(f"    η_μ/η_e = {eta_mu/eta_e:.8f}")
        print(f"    η_τ/η_e = {eta_tau/eta_e:.8f}")

    except Exception as e:
        print(f"    Failed to find g₀^τ: {e}")

print(f"\n{'='*78}")
