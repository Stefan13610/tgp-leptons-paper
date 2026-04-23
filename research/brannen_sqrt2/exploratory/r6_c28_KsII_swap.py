#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c28_KsII_swap.py:
Compute K_s^(II) via Fubini swap, analogous to K_c^(II) in r6_c6.

K_s^(II) = int_0^inf sin(s)*s*f'(s)*h'(s) ds = -B1-B2+B3-B4

where B_i = int_0^inf tau_i(t)*A(t)*Psi_i(t) dt
with A(t) = t*(f'(t))^2, and:

tau_1 = cos(t), Psi_1(t) = int_t^inf sin(s)cos(s)*f'(s) ds = Phi_2(t)  [reuse!]
tau_2 = sin(t), Psi_2(t) = int_t^inf sin^2(s)*f'(s) ds = -sin(t)/t - Phi_1(t)
tau_3 = cos(t), Psi_3(t) = int_t^inf sin^2(s)*f'(s)/s ds
tau_4 = sin(t), Psi_4(t) = int_t^inf sin(s)cos(s)*f'(s)/s ds = Phi_3(t)  [reuse!]

Key relations:
  Phi_1+Psi_2 = int_t^inf f'(s)ds = -f(t) = -sin(t)/t
  Phi_4+Psi_3 = int_t^inf f'(s)/s ds = cos(t)/(2t) - sin(t)/(2t^2) - pi/4 + Si(t)/2

Closed forms:
  Psi_1(t) = Phi_2(t) = (Si(3t)-Si(t))/2 - (cos(t)-cos(3t))/(4t)
  Psi_2(t) = (-3sin(t)+sin(3t))/(4t) + (Ci(t)-Ci(3t))/2
  Psi_3(t) = (cos(3t)-cos(t))/(8t) + (-3sin(t)+sin(3t))/(8t^2) - pi/8 - Si(t)/8 + 3*Si(3t)/8
  Psi_4(t) = Phi_3(t) = (3sin(t)-sin(3t))/(8t) + 3(Ci(3t)-Ci(t))/8 + (cos(3t)-cos(t))/(8t^2)

Then P_sin = K_s^(I) - 2*K_s^(II) = pi/12 - 2*K_s^(II)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 40
pi = mp.pi
Ci = mp.ci
Si = mp.si

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*mp.cos(s) - mp.sin(s))/s**2

def A_kernel(t):
    return t * fp(t)**2

# ---- Phi functions from K_c^(II) (r6_c6) ----
def Phi1(t):
    return mp.mpf('-0.5')*Ci(t) + mp.mpf('0.5')*Ci(3*t) - (mp.sin(t)+mp.sin(3*t))/(4*t)

def Phi2(t):
    return (Si(3*t)-Si(t))/2 - (mp.cos(t)-mp.cos(3*t))/(4*t)

def Phi3(t):
    return (3*mp.sin(t)-mp.sin(3*t))/(8*t) + mp.mpf(3)/8*(Ci(3*t)-Ci(t)) + (mp.cos(3*t)-mp.cos(t))/(8*t**2)

def Phi4(t):
    return (5*mp.cos(t)-mp.cos(3*t))/(8*t) - pi/8 + mp.mpf(5)/8*Si(t) - mp.mpf(3)/8*Si(3*t) - (mp.sin(t)+mp.sin(3*t))/(8*t**2)

# ---- Psi functions for K_s^(II) ----
def Psi1(t):
    """= Phi_2(t)"""
    return Phi2(t)

def Psi2(t):
    """= -sin(t)/t - Phi_1(t)"""
    return -mp.sin(t)/t - Phi1(t)

def Psi3(t):
    """From Phi4+Psi3 = cos(t)/(2t) - sin(t)/(2t^2) - pi/4 + Si(t)/2"""
    return (mp.cos(3*t)-mp.cos(t))/(8*t) + (-3*mp.sin(t)+mp.sin(3*t))/(8*t**2) \
           - pi/8 - Si(t)/8 + 3*Si(3*t)/8

def Psi4(t):
    """= Phi_3(t)"""
    return Phi3(t)

print("="*78)
print("  K_s^(II) via Fubini swap (dps=40)")
print("="*78)

# Verify Psi functions numerically
print("\n  --- Verifying Psi functions at test points ---")
for t_val in [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(5)]:
    # Psi_2 numerical
    psi2_num = mp.quadosc(lambda s: mp.sin(s)**2 * fp(s), [t_val, mp.inf], period=2*pi)
    psi2_closed = Psi2(t_val)
    print(f"  t={mp.nstr(t_val,2)}: Psi2 num={mp.nstr(psi2_num,20)}, closed={mp.nstr(psi2_closed,20)}, diff={mp.nstr(psi2_num-psi2_closed,6)}")

    # Psi_3 numerical
    psi3_num = mp.quadosc(lambda s: mp.sin(s)**2 * fp(s) / s, [t_val, mp.inf], period=2*pi)
    psi3_closed = Psi3(t_val)
    print(f"       Psi3 num={mp.nstr(psi3_num,20)}, closed={mp.nstr(psi3_closed,20)}, diff={mp.nstr(psi3_num-psi3_closed,6)}")

# Compute B_i integrals
print("\n  --- Computing B_i = int tau_i(t)*A(t)*Psi_i(t) dt ---")

B1 = mp.quadosc(lambda t: mp.cos(t)*A_kernel(t)*Psi1(t), [0, mp.inf], omega=1)
print(f"  B1 = {mp.nstr(B1, 30)}")

B2 = mp.quadosc(lambda t: mp.sin(t)*A_kernel(t)*Psi2(t), [0, mp.inf], omega=1)
print(f"  B2 = {mp.nstr(B2, 30)}")

B3 = mp.quadosc(lambda t: mp.cos(t)*A_kernel(t)*Psi3(t), [0, mp.inf], omega=1)
print(f"  B3 = {mp.nstr(B3, 30)}")

B4 = mp.quadosc(lambda t: mp.sin(t)*A_kernel(t)*Psi4(t), [0, mp.inf], omega=1)
print(f"  B4 = {mp.nstr(B4, 30)}")

KsII = -B1 - B2 + B3 - B4
print(f"\n  K_s^(II) = -B1-B2+B3-B4 = {mp.nstr(KsII, 30)}")

# P_sin = K_s^(I) - 2*K_s^(II) = pi/12 - 2*K_s^(II)
KsI = pi/12
Psin = KsI - 2*KsII
print(f"\n  K_s^(I) = pi/12 = {mp.nstr(KsI, 30)}")
print(f"  P_sin = K_s^(I) - 2*K_s^(II) = {mp.nstr(Psin, 30)}")

# Cross-check: verify K_c^(II) using the same framework
print(f"\n  --- Cross-check: K_c^(II) from r6_c6 ---")
A1 = mp.quadosc(lambda t: mp.cos(t)*A_kernel(t)*Phi1(t), [0, mp.inf], omega=1)
A2 = mp.quadosc(lambda t: mp.sin(t)*A_kernel(t)*Phi2(t), [0, mp.inf], omega=1)
A3 = mp.quadosc(lambda t: mp.cos(t)*A_kernel(t)*Phi3(t), [0, mp.inf], omega=1)
A4 = mp.quadosc(lambda t: mp.sin(t)*A_kernel(t)*Phi4(t), [0, mp.inf], omega=1)
KcII = -A1 - A2 + A3 - A4
KcII_known = mp.mpf('-0.031879037931728580134156')
print(f"  K_c^(II) = {mp.nstr(KcII, 25)}")
print(f"  known    = {mp.nstr(KcII_known, 25)}")
print(f"  diff     = {mp.nstr(KcII - KcII_known, 8)}")

# Now compute the perturbation amplitudes
print(f"\n{'='*78}")
print(f"  TAIL AMPLITUDES for alpha_5 formula")
print(f"{'='*78}")

ln2 = mp.log(2)
ln3 = mp.log(3)

B2_val = mp.mpf(1)/2 - ln3/8  # I_cos (proven)
A2_val = pi/8                  # pi/8 (proven)
B3_val = mp.mpf('0.012615939290114711837850')  # P_cos (30 digits)
A3_val = -Psin  # cos-amplitude of u3 (NEW!)

print(f"  B2 (I_cos, proven)   = {mp.nstr(B2_val, 25)}")
print(f"  A2 (pi/8, proven)    = {mp.nstr(A2_val, 25)}")
print(f"  B3 (P_cos, 30 dig)   = {mp.nstr(B3_val, 25)}")
print(f"  A3 = -P_sin (NEW!)   = {mp.nstr(A3_val, 25)}")

# Verify alpha_3 = B3 + A2^2/2
alpha3_check = B3_val + A2_val**2/2
alpha3_known = mp.mpf('0.08972222367362532604749')
print(f"\n  Verify: alpha_3 = B3 + A2^2/2 = {mp.nstr(alpha3_check, 25)}")
print(f"  known: {mp.nstr(alpha3_known, 25)}")
print(f"  diff:  {mp.nstr(alpha3_check - alpha3_known, 8)}")

# Partial alpha_5 (without A4, B4 which need f4 equation)
# alpha_5 = A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
# Terms not involving A4: A3^2/2 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
a5_partial = A3_val**2/2 - B2_val*A2_val*A3_val + A2_val**2*(B2_val**2-B3_val)/2 - A2_val**4/8
print(f"\n  Partial alpha_5 (without A2*A4 term):")
print(f"    = A3^2/2 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8")
print(f"    = {mp.nstr(a5_partial, 20)}")
print(f"  Total alpha_5 = {mp.nstr(a5_partial, 10)} + A2*A4")
print(f"  A2*A4 needed; from ODE: alpha_5 ≈ 0.0275")
print(f"  So A2*A4 ≈ 0.0275 - ({float(a5_partial):.6f}) = {0.0275 - float(a5_partial):.6f}")
print(f"  → A4 ≈ {(0.0275 - float(a5_partial))/float(A2_val):.6f}")

print("="*78)
