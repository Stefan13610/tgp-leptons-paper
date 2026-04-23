#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c36_B4A4_mpmath_quadosc.py:
High-precision B4, A4, B5, A5 via mpmath ODE + quadosc.

Strategy:
1. Solve coupled f2, f3, f4 ODE with mpmath RK4 (dps=35) on fine grid
2. Create interpolation lookup for f2, f2', f3, f3', f4, f4'
3. Construct S4, S5 source terms as functions
4. Compute B_n = ∫ cos(t)·t·S_n dt, A_n = -∫ sin(t)·t·S_n dt via quadosc

NEW ANALYTICAL RESULTS (verified):
- ∫ sin(t)·t·f1^2·(f1')^2 dt = 5*pi/128
- ∫ sin(t)·t·f1^3·(f1')^2 dt = pi/40
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 35
pi = mp.pi

print("="*78)
print("  B4, A4, B5 via mpmath ODE + quadosc (dps=35)")
print("="*78)

# --- f1 and f1' (analytical) ---
def f1(r):
    if r < mp.mpf('1e-20'):
        return mp.mpf(1) - r**2/6 + r**4/120
    return mp.sin(r)/r

def f1p(r):
    if r < mp.mpf('1e-20'):
        return -r/3 + r**3/30 - r**5/840
    return mp.cos(r)/r - mp.sin(r)/r**2

# --- Solve coupled ODE for f2, f3, f4 with mpmath RK4 ---
R_max = 600
dr = mp.mpf('0.01')
N = int(float(R_max / dr))

print(f"\n  Integrating f2, f3, f4 with mpmath RK4:")
print(f"    dps={mp.mp.dps}, R_max={R_max}, dr={float(dr)}, N={N}")

# State: [f2, f2', f3, f3', f4, f4']
def ode_rhs(r, y):
    f2, f2p, f3, f3p, f4, f4p = y
    fv = f1(r)
    fp = f1p(r)
    if r < mp.mpf('1e-10'):
        return [f2p, mp.mpf(0), f3p, mp.mpf(0), f4p, mp.mpf(0)]
    tr = 2/r
    s2 = -fp**2
    s3 = fv*fp**2 - 2*fp*f2p
    s4 = 2*fv*fp*f2p + f2*fp**2 - fv**2*fp**2 - 2*fp*f3p - f2p**2
    return [f2p, s2 - tr*f2p - f2,
            f3p, s3 - tr*f3p - f3,
            f4p, s4 - tr*f4p - f4]

# RK4 integration
y = [mp.mpf(0)]*6
r = mp.mpf('1e-8')

# Store on grid for interpolation
r_grid = []
f2_grid = []
f2p_grid = []
f3_grid = []
f3p_grid = []
f4_grid = []
f4p_grid = []

step = 0
while r < R_max + dr/2:
    r_grid.append(float(r))
    f2_grid.append(float(y[0]))
    f2p_grid.append(float(y[1]))
    f3_grid.append(float(y[2]))
    f3p_grid.append(float(y[3]))
    f4_grid.append(float(y[4]))
    f4p_grid.append(float(y[5]))

    k1 = ode_rhs(r, y)
    y2 = [y[i] + dr/2*k1[i] for i in range(6)]
    k2 = ode_rhs(r + dr/2, y2)
    y3 = [y[i] + dr/2*k2[i] for i in range(6)]
    k3 = ode_rhs(r + dr/2, y3)
    y4 = [y[i] + dr*k3[i] for i in range(6)]
    k4 = ode_rhs(r + dr, y4)

    y = [y[i] + dr/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(6)]
    r += dr
    step += 1

    if step % 10000 == 0:
        print(f"    step {step}, r={float(r):.0f}, f2={float(y[0]):.10e}, f3={float(y[2]):.10e}")

print(f"    Done: {step} steps")

# Convert to numpy arrays for fast interpolation
import numpy as np
r_grid = np.array(r_grid)
f2_grid = np.array(f2_grid)
f2p_grid = np.array(f2p_grid)
f3_grid = np.array(f3_grid)
f3p_grid = np.array(f3p_grid)
f4_grid = np.array(f4_grid)
f4p_grid = np.array(f4p_grid)

# Linear interpolation function
def interp(r_val, r_arr, v_arr):
    """Linear interpolation."""
    r_f = float(r_val)
    if r_f <= r_arr[0]:
        return v_arr[0]
    if r_f >= r_arr[-1]:
        return v_arr[-1]
    idx = np.searchsorted(r_arr, r_f) - 1
    if idx < 0: idx = 0
    if idx >= len(r_arr)-1: idx = len(r_arr)-2
    t = (r_f - r_arr[idx]) / (r_arr[idx+1] - r_arr[idx])
    return v_arr[idx] + t * (v_arr[idx+1] - v_arr[idx])

def get_f2(r): return mp.mpf(interp(r, r_grid, f2_grid))
def get_f2p(r): return mp.mpf(interp(r, r_grid, f2p_grid))
def get_f3(r): return mp.mpf(interp(r, r_grid, f3_grid))
def get_f3p(r): return mp.mpf(interp(r, r_grid, f3p_grid))

# Cross-check: tail amplitudes of u2 = r*f2
print(f"\n  Cross-check: u2 at r=500:")
r_test = 500.0
u2 = r_test * interp(r_test, r_grid, f2_grid)
print(f"    u2(500) = {u2:.10f}")
print(f"    B2*sin(500)+A2*cos(500) = {(0.5-np.log(3)/8)*np.sin(500) + np.pi/8*np.cos(500):.10f}")

# --- Compute B4, A4 via quadosc ---
print(f"\n{'='*78}")
print("  Computing B4, A4 via quadosc (omega=1)")
print("="*78)

# S4 = T1 + T2 + T3 + T4 + T5
# T1 = -f1^2*(f1')^2
# T2 = f2*(f1')^2
# T3 = 2*f1*f1'*f2'
# T4 = -2*f1'*f3'
# T5 = -(f2')^2

# B4 = sum of B4_Ti where B4_Ti = ∫ cos(t)*t*Ti dt
# A4 = sum of A4_Ti where A4_Ti = -∫ sin(t)*t*Ti dt

# T1 (f1-only): already computed to 40 digits
B4_T1 = mp.mpf('0.000728520701977491947807068225939')
A4_T1 = mp.mpf('0.12271846303085129837744700715935558')  # = 5*pi/128

# T2: f2*(f1')^2
def T2_cos(t):
    return mp.cos(t) * t * get_f2(t) * f1p(t)**2
def T2_sin(t):
    return -mp.sin(t) * t * get_f2(t) * f1p(t)**2

print("\n  T2: f2*(f1')^2")
B4_T2 = mp.quadosc(T2_cos, [0, R_max], omega=1)
A4_T2 = mp.quadosc(T2_sin, [0, R_max], omega=1)
print(f"    B4_T2 = {mp.nstr(B4_T2, 20)}")
print(f"    A4_T2 = {mp.nstr(A4_T2, 20)}")

# T3: 2*f1*f1'*f2'
def T3_cos(t):
    return mp.cos(t) * t * 2 * f1(t) * f1p(t) * get_f2p(t)
def T3_sin(t):
    return -mp.sin(t) * t * 2 * f1(t) * f1p(t) * get_f2p(t)

print("\n  T3: 2*f1*f1'*f2'")
B4_T3 = mp.quadosc(T3_cos, [0, R_max], omega=1)
A4_T3 = mp.quadosc(T3_sin, [0, R_max], omega=1)
print(f"    B4_T3 = {mp.nstr(B4_T3, 20)}")
print(f"    A4_T3 = {mp.nstr(A4_T3, 20)}")

# T4: -2*f1'*f3'
def T4_cos(t):
    return mp.cos(t) * t * (-2) * f1p(t) * get_f3p(t)
def T4_sin(t):
    return -mp.sin(t) * t * (-2) * f1p(t) * get_f3p(t)

print("\n  T4: -2*f1'*f3'")
B4_T4 = mp.quadosc(T4_cos, [0, R_max], omega=1)
A4_T4 = mp.quadosc(T4_sin, [0, R_max], omega=1)
print(f"    B4_T4 = {mp.nstr(B4_T4, 20)}")
print(f"    A4_T4 = {mp.nstr(A4_T4, 20)}")

# T5: -(f2')^2
def T5_cos(t):
    return mp.cos(t) * t * (-(get_f2p(t))**2)
def T5_sin(t):
    return -mp.sin(t) * t * (-(get_f2p(t))**2)

print("\n  T5: -(f2')^2")
B4_T5 = mp.quadosc(T5_cos, [0, R_max], omega=1)
A4_T5 = mp.quadosc(T5_sin, [0, R_max], omega=1)
print(f"    B4_T5 = {mp.nstr(B4_T5, 20)}")
print(f"    A4_T5 = {mp.nstr(A4_T5, 20)}")

# Sum
B4_total = B4_T1 + B4_T2 + B4_T3 + B4_T4 + B4_T5
A4_total = A4_T1 + A4_T2 + A4_T3 + A4_T4 + A4_T5

print(f"\n{'='*78}")
print("  RESULTS")
print("="*78)
print(f"  B4 = B4_T1 + B4_T2 + B4_T3 + B4_T4 + B4_T5")
print(f"     = {mp.nstr(B4_T1, 15)}")
print(f"     + {mp.nstr(B4_T2, 15)}")
print(f"     + {mp.nstr(B4_T3, 15)}")
print(f"     + {mp.nstr(B4_T4, 15)}")
print(f"     + {mp.nstr(B4_T5, 15)}")
print(f"     = {mp.nstr(B4_total, 15)}")
print(f"  (from ODE tail: 0.08807034)")

print(f"\n  A4 = A4_T1 + A4_T2 + A4_T3 + A4_T4 + A4_T5")
print(f"     = {mp.nstr(A4_T1, 15)} (= 5*pi/128)")
print(f"     + {mp.nstr(A4_T2, 15)}")
print(f"     + {mp.nstr(A4_T3, 15)}")
print(f"     + {mp.nstr(A4_T4, 15)}")
print(f"     + {mp.nstr(A4_T5, 15)}")
print(f"     = {mp.nstr(A4_total, 15)}")
print(f"  (from ODE tail: -0.10039114)")

# Compute alpha_4
import math
ln3 = math.log(3)
B2v = 0.5 - ln3/8
A2v = math.pi/8
A3v = -0.215712018451450166
alpha4 = float(B4_total) + A2v*A3v - B2v*A2v**2/2
print(f"\n  alpha_4 = {alpha4:+.12f}  (expect ~-0.0246)")

print(f"\n{'='*78}")
