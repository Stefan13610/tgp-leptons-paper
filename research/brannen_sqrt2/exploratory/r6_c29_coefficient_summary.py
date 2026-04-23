#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c29_coefficient_summary.py:
Summary of all known perturbation coefficients and analytical structures.

Full expansion: eta(delta) = 1 + alpha_2*d + alpha_3*d^2 + alpha_4*d^3 + alpha_5*d^4 + ...

Known:
  alpha_2 = c1/2 = (1 - ln3/4)/2 = I_cos (PROVEN)
  alpha_3 = pi^2/128 + P_cos = 0.089722... (30 digits, no closed form for P_cos)
  alpha_4 ≈ -0.0246 (from ODE, not yet computed via perturbation)
  alpha_5 ≈ 0.0275 (from ODE, 4 digits)

Tail amplitudes (u_n = r*f_n ~ B_n*sin(r) + A_n*cos(r)):
  B1 = 1, A1 = 0
  B2 = I_cos = 1/2 - ln3/8 (PROVEN)
  A2 = pi/8 (PROVEN)
  B3 = P_cos = 0.012616... (30 digits)
  A3 = -P_sin = -0.215712... (25 digits, NEW)

Analytical structures:
  K_c^(I) = (ln2-1)/6 (PROVEN)
  K_s^(I) = pi/12 (PROVEN)
  K_c^(II) = -0.031879... (30 digits, no closed form)
  K_s^(II) = 0.023044... (25 digits, no closed form)

This script:
1. Verifies all known values
2. Computes alpha_4 from perturbation amplitudes
3. Tests PSLQ on P_sin, K_s^(II)
4. Checks for patterns between the "new" constants
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 35
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

print("="*78)
print("  PERTURBATION COEFFICIENT SUMMARY")
print("="*78)

# Known amplitudes
B1 = mp.mpf(1)
A1 = mp.mpf(0)
B2 = mp.mpf(1)/2 - ln3/8
A2 = pi/8
B3 = mp.mpf('0.012615939290114711837850')
A3 = mp.mpf('-0.215712018451450166339614')  # -P_sin
B4_ode = mp.mpf('0.088')  # rough from ODE, not reliable
A4_ode = mp.mpf('-0.098')  # rough from ODE, not reliable

# Derived coefficients
alpha2 = B2  # = I_cos
alpha3 = B3 + A2**2/2  # = P_cos + pi^2/128
alpha4_pert = B4_ode + A2*A3 - B2*A2**2/2  # rough, from formula
# alpha_5 requires B5 + other terms

print(f"\n  Tail amplitudes:")
print(f"    B1 = {mp.nstr(B1, 20)}")
print(f"    A1 = {mp.nstr(A1, 20)}")
print(f"    B2 = I_cos = {mp.nstr(B2, 25)} (PROVEN)")
print(f"    A2 = pi/8 = {mp.nstr(A2, 25)} (PROVEN)")
print(f"    B3 = P_cos = {mp.nstr(B3, 25)} (30 digits)")
print(f"    A3 = -P_sin = {mp.nstr(A3, 25)} (25 digits, NEW)")

print(f"\n  Eta expansion coefficients:")
print(f"    alpha_2 = c1/2 = {mp.nstr(alpha2, 20)} (PROVEN)")
print(f"    alpha_3 = {mp.nstr(alpha3, 20)} (30 digits)")
print(f"    alpha_4 ~ {mp.nstr(alpha4_pert, 8)} (rough from ODE)")
print(f"    alpha_5 ~ 0.0275 (4 digits from ODE)")

# Known analytical constants
KcI = (ln2 - 1) / 6
KsI = pi / 12
KcII = mp.mpf('-0.031879037931728580134156')
KsII = mp.mpf('0.0230436846738496350994699')

print(f"\n  K-constants:")
print(f"    K_c^(I) = (ln2-1)/6 = {mp.nstr(KcI, 25)} (PROVEN)")
print(f"    K_s^(I) = pi/12 = {mp.nstr(KsI, 25)} (PROVEN)")
print(f"    K_c^(II) = {mp.nstr(KcII, 25)} (30 digits)")
print(f"    K_s^(II) = {mp.nstr(KsII, 25)} (25 digits)")

# P_cos, P_sin
Pcos = KcI - 2*KcII
Psin = KsI - 2*KsII
print(f"\n  P-constants:")
print(f"    P_cos = K_c^(I) - 2*K_c^(II) = {mp.nstr(Pcos, 25)}")
print(f"    P_sin = K_s^(I) - 2*K_s^(II) = {mp.nstr(Psin, 25)}")
print(f"    P_cos/P_sin = {mp.nstr(Pcos/Psin, 15)}")
print(f"    arctan(P_sin/P_cos) = {mp.nstr(mp.atan2(Psin, Pcos), 15)} rad")
print(f"    = {mp.nstr(mp.atan2(Psin, Pcos)*180/pi, 10)} deg")

# PSLQ tests on P_sin and K_s^(II)
print(f"\n{'='*78}")
print("  PSLQ TESTS on NEW constants")
print("="*78)

# Test P_sin
targets = {
    'P_sin': Psin,
    'K_s^(II)': KsII,
    'K_c^(II)': KcII,
}

basis_names = ['pi^2', 'pi^4', 'ln2', 'ln3', 'ln2^2', 'ln3^2', 'ln2*ln3',
               'zeta3', 'G', 'pi*ln2', 'pi*ln3', '1']
basis_vals = [pi**2, pi**4, ln2, ln3, ln2**2, ln3**2, ln2*ln3,
              mp.zeta(3), mp.catalan, pi*ln2, pi*ln3, mp.mpf(1)]

for tname, tval in targets.items():
    print(f"\n  --- PSLQ for {tname} = {mp.nstr(tval, 20)} ---")
    basis = [tval] + basis_vals
    names = [tname] + basis_names
    try:
        rel = mp.pslq(basis, maxcoeff=10**8, tol=mp.mpf(10)**(-20))
        if rel and rel[0] != 0:
            parts = []
            for n, c in zip(names[1:], rel[1:]):
                if c != 0:
                    parts.append(f"{c:+d}*{n}")
            print(f"    {rel[0]}*{tname} = {', '.join(parts)}")
            check = sum(c*v for c,v in zip(rel, basis))
            print(f"    residual = {mp.nstr(check, 6)}")
        elif rel:
            print(f"    trivial (coeff of target = 0)")
        else:
            print(f"    no relation found")
    except Exception as e:
        print(f"    error: {e}")

# Test relationships BETWEEN the unknown constants
print(f"\n  --- Relations between unknown constants ---")

# Are K_c^(II) and K_s^(II) related?
rel_tests = [
    ('K_c^(II)/K_s^(II)', KcII/KsII),
    ('K_c^(II) + K_s^(II)', KcII + KsII),
    ('K_c^(II) - K_s^(II)', KcII - KsII),
    ('K_c^(II)^2 + K_s^(II)^2', KcII**2 + KsII**2),
    ('P_cos/P_sin', Pcos/Psin),
    ('P_cos^2 + P_sin^2', Pcos**2 + Psin**2),
]

for name, val in rel_tests:
    ident = mp.identify(val)
    print(f"  {name:30s} = {mp.nstr(val, 15):>20s}  identify: {ident}")

# Test: is P_cos^2 + P_sin^2 a nice constant?
Pmod2 = Pcos**2 + Psin**2
print(f"\n  P_cos^2 + P_sin^2 = {mp.nstr(Pmod2, 20)}")
basis_p = [Pmod2, pi**2, pi**4, ln2, ln3, ln2**2, mp.zeta(3), mp.mpf(1)]
names_p = ['|P|^2', 'pi^2', 'pi^4', 'ln2', 'ln3', 'ln2^2', 'zeta3', '1']
try:
    rel = mp.pslq(basis_p, maxcoeff=10**8)
    if rel and rel[0] != 0:
        parts = [f"{c:+d}*{n}" for n,c in zip(names_p[1:], rel[1:]) if c != 0]
        print(f"  PSLQ: {rel[0]}*|P|^2 = {', '.join(parts)}")
    else:
        print(f"  PSLQ: no relation for |P|^2")
except:
    print(f"  PSLQ failed")

print("="*78)
