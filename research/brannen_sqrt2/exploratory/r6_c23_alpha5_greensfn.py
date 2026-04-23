#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c23_alpha5_greensfn.py:
High-precision alpha_5 via Green's function integral accumulation.

Key insight: instead of fitting oscillatory tails (inaccurate at finite R),
accumulate the Green's function integrals directly:
  B_n = int_0^inf cos(s) * source_n(s) ds  (sin-amplitude of u_n)
  A_n = -int_0^inf sin(s) * source_n(s) ds  (cos-amplitude)

These converge monotonically (no oscillation fitting needed).

We integrate the ODE system AND the tail integrals simultaneously.
The 4th-order formula (derived in r6_c22):
  alpha_5 = q^2/2 + p*r - a*p*q + p^2*(a^2-b)/2 - p^4/8
where a=B2, b=B3, c=B4, p=A2, q=A3, r=A4.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 30
pi = mp.pi

print("="*78)
print("  ALPHA_5 via Green's function integral accumulation")
print("  dps =", mp.mp.dps)
print("="*78)

def f1(r):
    if r < mp.mpf('1e-30'): return mp.mpf(1)
    return mp.sin(r) / r

def f1p(r):
    if r < mp.mpf('1e-30'): return mp.mpf(0)
    return mp.cos(r)/r - mp.sin(r)/r**2

R_MAX = 400
N_steps = 160000  # dr = 0.0025
dr = mp.mpf(R_MAX) / N_steps

print(f"\n  R_max = {R_MAX}, N_steps = {N_steps}, dr = {mp.nstr(dr, 6)}")

# State vector: [u2, u2', u3, u3', u4, u4',
#                Bc2, Bs2,   (accum integrals for u2 tail)
#                Bc3, Bs3,   (accum integrals for u3 tail)
#                Bc4, Bs4]   (accum integrals for u4 tail)
# where Bc_n = int_0^r cos(s)*src_n(s) ds
#       Bs_n = int_0^r sin(s)*src_n(s) ds

def derivs(r, st):
    u2, u2p, u3, u3p, u4, u4p = st[:6]

    f1v = f1(r)
    f1pv = f1p(r)

    # Source for u2: -r*(f1')^2
    src2 = -r * f1pv**2

    # f2, f2' from u2
    if r < mp.mpf('1e-15'):
        f2v = mp.mpf(0)
        f2pv = mp.mpf(0)
    else:
        f2v = u2 / r
        f2pv = (u2p * r - u2) / r**2

    # Source for u3: r*[-(2*f1'*f2' - f1*(f1')^2)]
    src3 = r * (-(2*f1pv*f2pv - f1v*f1pv**2))

    # f3, f3' from u3
    if r < mp.mpf('1e-15'):
        f3v = mp.mpf(0)
        f3pv = mp.mpf(0)
    else:
        f3v = u3 / r
        f3pv = (u3p * r - u3) / r**2

    # Source for u4: r*[-((f2')^2 + 2*f1'*f3' - 2*f1*f1'*f2' - f2*(f1')^2 + f1^2*(f1')^2)]
    src4 = r * (-((f2pv)**2 + 2*f1pv*f3pv - 2*f1v*f1pv*f2pv - f2v*f1pv**2 + f1v**2*f1pv**2))

    cosr = mp.cos(r)
    sinr = mp.sin(r)

    return [
        u2p, -u2 + src2,
        u3p, -u3 + src3,
        u4p, -u4 + src4,
        cosr * src2, sinr * src2,   # d/dr of Bc2, Bs2
        cosr * src3, sinr * src3,   # d/dr of Bc3, Bs3
        cosr * src4, sinr * src4,   # d/dr of Bc4, Bs4
    ]

# RK4 stepper
def rk4_step(r, state, h):
    k1 = derivs(r, state)
    s1 = [state[i] + h/2*k1[i] for i in range(12)]
    k2 = derivs(r + h/2, s1)
    s2 = [state[i] + h/2*k2[i] for i in range(12)]
    k3 = derivs(r + h/2, s2)
    s3 = [state[i] + h*k3[i] for i in range(12)]
    k4 = derivs(r + h, s3)
    return [state[i] + h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(12)]

# Initialize
r = mp.mpf('1e-8')
state = [mp.mpf(0)] * 12  # all zeros

print("  Integrating...")
checkpoints = [50, 100, 150, 200, 300, 400]
check_idx = 0

for step in range(N_steps):
    state = rk4_step(r, state, dr)
    r += dr

    if check_idx < len(checkpoints) and float(r) >= checkpoints[check_idx]:
        Bc2, Bs2 = state[6], state[7]
        Bc3, Bs3 = state[8], state[9]
        Bc4, Bs4 = state[10], state[11]

        B2 = Bc2   # sin-amp of u2
        A2 = -Bs2  # cos-amp of u2
        B3 = Bc3
        A3 = -Bs3
        B4 = Bc4
        A4 = -Bs4

        # Compute alpha_3 and alpha_5 at this checkpoint
        a3_here = B3 + A2**2/2
        a5_here = A3**2/2 + A2*A4 - B2*A2*A3 + A2**2*(B2**2 - B3)/2 - A2**4/8

        print(f"\n    r = {checkpoints[check_idx]}:")
        print(f"      B2 = {mp.nstr(B2, 20)}  (I_cos, expected 0.36267346...)")
        print(f"      A2 = {mp.nstr(A2, 20)}  (pi/8, expected 0.39269908...)")
        print(f"      B3 = {mp.nstr(B3, 20)}  (P_cos, expected 0.01261594...)")
        print(f"      A3 = {mp.nstr(A3, 20)}  (-P_sin)")
        print(f"      B4 = {mp.nstr(B4, 20)}  (Q_cos)")
        print(f"      A4 = {mp.nstr(A4, 20)}  (-Q_sin)")
        print(f"      alpha_3 = {mp.nstr(a3_here, 15)}  (expected 0.089722224...)")
        print(f"      alpha_5 = {mp.nstr(a5_here, 15)}")

        check_idx += 1

# Final values
Bc2, Bs2 = state[6], state[7]
Bc3, Bs3 = state[8], state[9]
Bc4, Bs4 = state[10], state[11]

B2_final = Bc2
A2_final = -Bs2
B3_final = Bc3
A3_final = -Bs3
B4_final = Bc4
A4_final = -Bs4

a = B2_final
b = B3_final
c = B4_final
p = A2_final
q = A3_final
rv = A4_final

alpha2 = a
alpha3 = b + p**2/2
alpha4 = c + p*q - a*p**2/2
alpha5 = q**2/2 + p*rv - a*p*q + p**2*(a**2 - b)/2 - p**4/8

print(f"\n{'='*78}")
print(f"  FINAL RESULTS (R_max = {R_MAX}):")
print(f"{'='*78}")
print(f"  Tail amplitudes:")
print(f"    B2 (I_cos)  = {mp.nstr(B2_final, 25)}")
print(f"    A2 (pi/8)   = {mp.nstr(A2_final, 25)}")
print(f"    B3 (P_cos)  = {mp.nstr(B3_final, 25)}")
print(f"    A3 (-P_sin) = {mp.nstr(A3_final, 25)}")
print(f"    B4 (Q_cos)  = {mp.nstr(B4_final, 25)}")
print(f"    A4 (-Q_sin) = {mp.nstr(A4_final, 25)}")

print(f"\n  Perturbation coefficients:")
print(f"    alpha_2 = {mp.nstr(alpha2, 20)}  (expected c1/2 = {mp.nstr(mp.mpf(1)/2 - mp.log(3)/8, 20)})")
print(f"    alpha_3 = {mp.nstr(alpha3, 20)}  (expected 0.089722223674...)")
print(f"    alpha_4 = {mp.nstr(alpha4, 20)}  (new!)")
print(f"    alpha_5 = {mp.nstr(alpha5, 20)}  (compare ODE: ~0.0275)")

# Known values for comparison
B2_exact = mp.mpf(1)/2 - mp.log(3)/8
A2_exact = pi/8
B3_exact = mp.mpf('0.012615939290114711837850')
alpha3_exact = mp.mpf('0.08972222367362532604749')

print(f"\n  Errors (from known 30-digit values):")
print(f"    B2 err = {mp.nstr(B2_final - B2_exact, 6)}")
print(f"    A2 err = {mp.nstr(A2_final - A2_exact, 6)}")
print(f"    B3 err = {mp.nstr(B3_final - B3_exact, 6)}")
print(f"    alpha_3 err = {mp.nstr(alpha3 - alpha3_exact, 6)}")

# Check K_s^(I) = pi/12 conjecture
KsI = mp.quadosc(lambda s: mp.sin(s) * s * f1(s) * f1p(s)**2, [0, mp.inf], omega=1)
print(f"\n  Bonus: K_s^(I) = {mp.nstr(KsI, 25)}")
print(f"         pi/12   = {mp.nstr(pi/12, 25)}")
print(f"         diff    = {mp.nstr(KsI - pi/12, 8)}")

print("="*78)
