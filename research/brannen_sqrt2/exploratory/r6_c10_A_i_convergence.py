#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c10_A_i_convergence.py: Investigate convergence of A_i integrands at t=0.

For the A_i = int cos/sin(t) * A(t) * Phi_i(t) dt, A(t) ~ t^3/9 at t=0.
So even if Phi_i(t) diverges like 1/t at t=0, A(t)·Phi_i(t) ~ t^2, integrable.

The issue is when we DECOMPOSE A_i into sums of integrals by expanding A(t)·Phi_i(t)
into elementary pieces — individual pieces may diverge, but the SUM converges.

Examples:
   cos^2(s)f'(s) = (3cos(s) + cos(3s))/(4s) - (sin(s) + sin(3s))/(4s^2)
   Near s=0: f'(s) ~ -s/3, so cos^2 f'(s) ~ -s/3. Finite. But the RHS pieces individually
   behave: (3+1)/(4s) = 1/s  (divergent) and -(s+3s)/(4s^2) = -1/s  (divergent).
   Sum: 1/s - 1/s = 0 at leading order. Next: from (1/4s)·(3(1-s^2/2)+(1-9s^2/2)) - (1/4s^2)·(4s - (s^3/6 + 27s^3/6))
       = (1/4s)·(4 - 6s^2) - (1/4s^2)·(4s - 14s^3/3) + O(s^5)
       = 1/s - 3s/2 - 1/s + 7s/6 = -3s/2 + 7s/6 = -4s/12 = -s/3. ✓

So when we decompose the A_i, individual pieces DIVERGE at t=0, requiring regularization
or careful combination. The ChiFunctional (chi_2) approach above is fine as long as we
keep things combined.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
Ci = mp.ci
Si = mp.si
ln2 = mp.log(2)
ln3 = mp.log(3)

def f_fun(s):
    if s < mp.mpf('1e-20'): return mp.mpf(1) - s**2/6
    return mp.sin(s)/s

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360
    return (s*mp.cos(s) - mp.sin(s))/s**2

def A_ker(t): return t * fp(t)**2

def Phi1(t):
    return -Ci(t)/2 + Ci(3*t)/2 - (mp.sin(t) + mp.sin(3*t))/(4*t)

# Test: A(t) * cos(t) * Phi_1(t) ~ t^2 * something at t=0.
print("  A(t)·cos(t)·Phi_1(t) at small t (should be ~ O(t^2)):")
for t in [0.001, 0.01, 0.1, 1.0]:
    t_mp = mp.mpf(t)
    val = A_ker(t_mp) * mp.cos(t_mp) * Phi1(t_mp)
    print(f"    t={t:5.3f}: A(t)*cos(t)*Phi_1(t) = {mp.nstr(val, 10):>18}  /t^2 = {mp.nstr(val/t_mp**2, 10)}")

# Alternative: cos(t) t(f')^2 expanded to frequency components.
# cos^2(t)/t near t=0 ~ 1/t (divergent).
# -sin(2t)/t^2 near t=0 ~ -2/t (divergent).
# +sin^2(t)/t^3 near t=0 ~ t^2/t^3 = 1/t.
# Sum: 1/t - 2/t + 1/t = 0. Good, cancellation.
# But individual pieces diverge.

print("\n  t(f')^2 vs frequency decomposition at small t:")
for t in [0.001, 0.01, 0.1, 1.0]:
    t_mp = mp.mpf(t)
    fp_v = fp(t_mp)
    A_direct = t_mp * fp_v**2

    # Decomposed:
    p1 = mp.cos(t_mp)**2 / t_mp
    p2 = -mp.sin(2*t_mp) / t_mp**2
    p3 = mp.sin(t_mp)**2 / t_mp**3
    A_decomp = p1 + p2 + p3

    print(f"    t={t:6.3f}: A(t)={mp.nstr(A_direct, 12):>18}  decomp={mp.nstr(A_decomp, 12):>18}  parts: {mp.nstr(p1,5)} + {mp.nstr(p2,5)} + {mp.nstr(p3,5)}")

# ==== Now apply frequency decomp of cos(t)·A(t) and Phi_1(t) ====
# cos(t)·A(t) = 3cos(t)/(4t) + cos(3t)/(4t) - [sin(t)+sin(3t)]/(2t^2) + [cos(t)-cos(3t)]/(4t^3)
# Phi_1(t)   = (1/2)[Ci(3t) - Ci(t)] - [sin(t)+sin(3t)]/(4t)

# Product: many cross terms. Let's index them:
#   P11 = (3cos(t)/(4t)) · (1/2)[Ci(3t)-Ci(t)] = 3cos(t)[Ci(3t)-Ci(t)]/(8t)
#   P12 = (cos(3t)/(4t)) · (1/2)[Ci(3t)-Ci(t)] = cos(3t)[Ci(3t)-Ci(t)]/(8t)
#   P13 = -(sin(t)+sin(3t))/(2t^2) · (1/2)[Ci(3t)-Ci(t)] = -(sin(t)+sin(3t))[Ci(3t)-Ci(t)]/(4t^2)
#   P14 = (cos(t)-cos(3t))/(4t^3) · (1/2)[Ci(3t)-Ci(t)] = (cos(t)-cos(3t))[Ci(3t)-Ci(t)]/(8t^3)
#   P21 = (3cos(t)/(4t)) · -(sin(t)+sin(3t))/(4t) = -3cos(t)(sin(t)+sin(3t))/(16t^2)
#   P22 = (cos(3t)/(4t)) · -(sin(t)+sin(3t))/(4t) = -cos(3t)(sin(t)+sin(3t))/(16t^2)
#   P23 = -(sin(t)+sin(3t))/(2t^2) · -(sin(t)+sin(3t))/(4t) = (sin(t)+sin(3t))^2/(8t^3)
#   P24 = (cos(t)-cos(3t))/(4t^3) · -(sin(t)+sin(3t))/(4t) = -(cos(t)-cos(3t))(sin(t)+sin(3t))/(16t^4)

# Each P_ij is a product of trigs · (1/t^n) · (possibly Ci(3t)-Ci(t)).
# Individual integrands P_ij(t) might diverge at 0, but their sum equals A(t)cos(t)Phi_1(t)
# which is O(t^2), integrable.

# Strategy: evaluate each P_ij via table integrals after further frequency decomposition.
# The key non-trivial integrals are those involving Ci(3t)-Ci(t) (the "log-Frullani" parts).

print(f"\n{'='*72}")
print(f"  A_1 = sum of 8 elementary integrals of products P_ij")
print(f"  All involve either pure trig/powers of t, or [Ci(3t)-Ci(t)]*trig/powers of t.")
print(f"  The latter are 'Frullani-type' and can be evaluated via derivatives under integral.")
print(f"{'='*72}")

# Pure trig/power integrals we need (reducing A_1 without Ci part):
# Sum P21 + P22 + P23 + P24 = cos(t)·A(t) · [-(sin(t)+sin(3t))/(4t)]
# But we can also express it directly since sin+sin / 4t is simple times cos(t)·A(t).
# Let me compute numerically to check decomposition.

# Define:
def P_pure(t):
    """cos(t)·A(t) · [-(sin(t)+sin(3t))/(4t)]  — the 'non-Ci' piece of Phi_1 times cos(t)·A(t)"""
    cos_tAt = mp.cos(t) * t * fp(t)**2
    return cos_tAt * (-(mp.sin(t) + mp.sin(3*t))/(4*t))

def P_Ci(t):
    """cos(t)·A(t) · (1/2)[Ci(3t)-Ci(t)]"""
    cos_tAt = mp.cos(t) * t * fp(t)**2
    return cos_tAt * (Ci(3*t) - Ci(t))/2

# Verify decomposition: P_pure(t) + P_Ci(t) = cos(t)A(t)Phi_1(t)
print(f"\n  Decomposition check P_pure + P_Ci = A_1 integrand:")
for t in [0.5, 1.0, 2.0, 5.0, 10.0]:
    t_mp = mp.mpf(t)
    full = mp.cos(t_mp) * A_ker(t_mp) * Phi1(t_mp)
    check = P_pure(t_mp) + P_Ci(t_mp)
    print(f"    t={t:5.2f}: full={mp.nstr(full, 12)} sum={mp.nstr(check, 12)} diff={mp.nstr(full-check, 3)}")

# Compute each piece
print(f"\n  Integrating A_1 in two pieces:")
A1_pure = mp.quadosc(P_pure, [0, mp.inf], period=2*pi)
A1_Ci = mp.quadosc(P_Ci, [0, mp.inf], period=2*pi)
A1_total = A1_pure + A1_Ci
print(f"    A_1 (pure, no Ci)  = {mp.nstr(A1_pure, 20)}")
print(f"    A_1 (Ci part)      = {mp.nstr(A1_Ci, 20)}")
print(f"    A_1 (sum)          = {mp.nstr(A1_total, 20)}")

# Expanded pure part:
# P_pure(t) = cos(t)·A(t) · [-(sin(t)+sin(3t))/(4t)]
# cos(t)A(t) = 3cos(t)/(4t) + cos(3t)/(4t) - (sin(t)+sin(3t))/(2t^2) + (cos(t)-cos(3t))/(4t^3)
# times -(sin(t)+sin(3t))/(4t):
#
# Term 1: -3cos(t)(sin(t)+sin(3t))/(16t^2)
#   = -3[sin(2t)/2 + (sin(4t)+sin(2t))/2]/(16t^2)    [cos(t)sin(t)=sin(2t)/2, cos(t)sin(3t)=(sin(4t)+sin(2t))/2]
#   = -3[sin(2t)/2 + sin(4t)/2 + sin(2t)/2]/(16t^2)
#   = -3[sin(2t) + sin(4t)/2]/(16t^2)
#   = -[3sin(2t)/16 + 3sin(4t)/32]/t^2
# Term 2: -cos(3t)(sin(t)+sin(3t))/(16t^2)
#   cos(3t)sin(t) = [sin(4t) - sin(2t)]/2, cos(3t)sin(3t) = sin(6t)/2
#   = -[sin(4t)-sin(2t) + sin(6t)]/(32t^2)
#   = [sin(2t) - sin(4t) - sin(6t)]/(32t^2)
# Term 3: +(sin(t)+sin(3t))^2/(8t^3)
#   (sin(t)+sin(3t))^2 = sin^2(t) + 2sin(t)sin(3t) + sin^2(3t)
#                      = (1-cos(2t))/2 + (cos(2t)-cos(4t)) + (1-cos(6t))/2
#                      = 1 + cos(2t)/2 - cos(4t) - cos(6t)/2
#   Hmm actually: (1-cos(2t))/2 + 2·(cos(2t)-cos(4t))/2 + (1-cos(6t))/2
#                = (1-cos(2t))/2 + (cos(2t)-cos(4t)) + (1-cos(6t))/2
#                = 1 + [-1/2 + 1]cos(2t) - cos(4t) - cos(6t)/2
#                = 1 + cos(2t)/2 - cos(4t) - cos(6t)/2
# Term 3 = [1 + cos(2t)/2 - cos(4t) - cos(6t)/2]/(8t^3)
# Term 4: -(cos(t)-cos(3t))(sin(t)+sin(3t))/(16t^4)
#   (cos(t)-cos(3t))(sin(t)+sin(3t)) = cos(t)sin(t) + cos(t)sin(3t) - cos(3t)sin(t) - cos(3t)sin(3t)
#     = sin(2t)/2 + [sin(4t)+sin(2t)]/2 - [sin(4t)-sin(2t)]/2 - sin(6t)/2
#     = sin(2t)/2 + sin(4t)/2 + sin(2t)/2 - sin(4t)/2 + sin(2t)/2 - sin(6t)/2
#     = 3sin(2t)/2 - sin(6t)/2
# Term 4 = -[3sin(2t)/2 - sin(6t)/2]/(16t^4) = -[3sin(2t) - sin(6t)]/(32t^4)

# So P_pure(t) = -3sin(2t)/(16t^2) - 3sin(4t)/(32t^2) + [sin(2t) - sin(4t) - sin(6t)]/(32t^2)
#              + [1 + cos(2t)/2 - cos(4t) - cos(6t)/2]/(8t^3)
#              - [3sin(2t) - sin(6t)]/(32t^4)
# Simplifying sin/t^2 coefficients:
#   sin(2t)/t^2: -3/16 + 1/32 = (-6+1)/32 = -5/32
#   sin(4t)/t^2: -3/32 - 1/32 = -4/32 = -1/8
#   sin(6t)/t^2: -1/32
# So sin/t^2 part: -5sin(2t)/(32t^2) - sin(4t)/(8t^2) - sin(6t)/(32t^2)
#
# 1/t^3 part: [1 + cos(2t)/2 - cos(4t) - cos(6t)/2]/(8t^3)
#
# 1/t^4 part: -[3sin(2t) - sin(6t)]/(32t^4)

# So P_pure has 4 integrals:
#  I_a = int -5sin(2t)/(32t^2) dt  -- individual diverges or is finite? sin(2t)/t^2 near 0 ~ 2/t diverges
#  Need to combine 1/t^2 + 1/t^3 + 1/t^4 terms together to get finite integrand at t=0.

# Verify P_pure is finite at t=0:
print(f"\n  P_pure(t) at small t (should be finite, ~ O(t^2) since A_1 integrand is O(t^3) and divided by t from Phi's 1/t):")
# Actually P_pure = cos(t)·A(t)·Phi_1_nonCi(t). cos(t)·A(t) ~ -t/3 at small t.
# Phi_1_nonCi(t) = -(sin(t)+sin(3t))/(4t) ~ -(t+3t)/(4t) = -1 at t=0.
# So P_pure(0+) = -t/3 · (-1) = t/3 → 0.
for t in [0.001, 0.01, 0.1, 1.0]:
    t_mp = mp.mpf(t)
    p = P_pure(t_mp)
    print(f"    t={t:6.3f}: P_pure(t) = {mp.nstr(p, 12)}")
print("="*72)
