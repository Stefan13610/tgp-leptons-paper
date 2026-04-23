#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c7b_KcII_direct.py: Compute K_c^(II) DIRECTLY via nested integrals (slow but reliable).
Purpose: verify whether the target (ln 2 - 1)/12 - 9 pi^2/14080 is exact or just close.

K_c^(II) = int_0^inf cos(s) s f'(s) h'(s) ds,  where:
    h'(s) = -u'(s)/s + u(s)/s^2
    u(s) = sin(s) J_c(s) - cos(s) J_s(s)
    u'(s) = cos(s) J_c(s) + sin(s) J_s(s)
    J_c(s) = int_0^s cos(t) t (f'(t))^2 dt
    J_s(s) = int_0^s sin(t) t (f'(t))^2 dt

Direct evaluation: pre-compute J_c, J_s on a dense grid, interpolate, integrate.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp, quad
mp.mp.dps = 30

pi_mp = mp.pi
ln2_mp = mp.log(2)

print("="*72)
print("  K_c^(II) direct nested evaluation")
print("="*72)

# Use ODE approach: dJ_c/dt = cos(t) t (f'(t))^2,  dJ_s/dt = sin(t) t (f'(t))^2
# Solve from t=0 to t=TMAX with tight tolerance.

def fp_np(s):
    small = s < 1e-3
    out_small = -s/3 + s**3/30 - s**5/840 + s**7/45360
    out_big = np.where(s > 0, (s*np.cos(s) - np.sin(s))/(s*s + 1e-300), out_small)
    return np.where(small, out_small, out_big)

def fp_scalar(s):
    if s < 1e-3:
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*np.cos(s) - np.sin(s))/(s*s)

def A_scalar(t):
    fpv = fp_scalar(t)
    return t * fpv * fpv

TMAX = 200.0  # Integrate 200 periods; oscillation damps (f')^2 ~ cos^2/s^2

def rhs(t, y):
    a = A_scalar(t)
    return [np.cos(t) * a, np.sin(t) * a]

print(f"  Solving ODE for J_c(t), J_s(t) on [0, {TMAX}]...")
sol = solve_ivp(rhs, (0.0, TMAX), [0.0, 0.0], method='DOP853', rtol=1e-13, atol=1e-15, dense_output=True)
print(f"    status: {sol.status}, msg: {sol.message}")
print(f"    J_c({TMAX}) = {sol.y[0, -1]:.15g}")
print(f"    J_s({TMAX}) = {sol.y[1, -1]:.15g}")

# Compare with theoretical J_c(inf) = ln(3)/8 - 1/2, J_s(inf) = pi/8
import math
Jc_inf = math.log(3)/8 - 0.5
Js_inf = math.pi/8
print(f"    J_c(inf) th = {Jc_inf:.15g}")
print(f"    J_s(inf) th = {Js_inf:.15g}")

# Build splines for u(s), u'(s)
t_dense = np.linspace(1e-6, TMAX, 200001)
yvals = sol.sol(t_dense)
Jc_arr = yvals[0]
Js_arr = yvals[1]

def u_arr(t, Jc, Js):
    return np.sin(t) * Jc - np.cos(t) * Js
def up_arr(t, Jc, Js):
    return np.cos(t) * Jc + np.sin(t) * Js

u_vals = u_arr(t_dense, Jc_arr, Js_arr)
up_vals = up_arr(t_dense, Jc_arr, Js_arr)

# h'(s) = -u'(s)/s + u(s)/s^2
hp_vals = -up_vals/t_dense + u_vals/t_dense**2

# Full integrand: cos(s) * s * f'(s) * h'(s)
fp_vals = fp_np(t_dense)
integrand = np.cos(t_dense) * t_dense * fp_vals * hp_vals

# Trapezoidal / Simpson integration
from scipy.integrate import simpson
KcII_direct = simpson(integrand, x=t_dense)
print(f"\n  K_c^(II) direct (simpson, TMAX={TMAX}) = {KcII_direct:.15g}")

# Target
KcII_target = float((ln2_mp - 1)/12 - 9*pi_mp**2/14080)
print(f"  K_c^(II) target (ln 2 - 1)/12 - 9 pi^2/14080 = {KcII_target:.15g}")
print(f"  diff = {KcII_direct - KcII_target:.5g}")

# Also try integrating [0, L] for several L and checking tail behavior
print(f"\n  Tail behavior check (cutoff L):")
for L in [50, 80, 100, 150, 200]:
    mask = t_dense <= L
    val = simpson(integrand[mask], x=t_dense[mask])
    print(f"    L={L:3d}: K_c^(II) = {val:.15g}  diff_target = {val - KcII_target:.5g}")

print("="*72)
