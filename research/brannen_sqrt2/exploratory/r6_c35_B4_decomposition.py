#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c35_B4_decomposition.py:
Decompose B4, A4 into structural components.

B4 = integral_0^inf cos(t)*t*S4(t) dt
A4 = -integral_0^inf sin(t)*t*S4(t) dt

where S4 = T1 + T2 + T3 + T4 + T5:
  T1 = -f1^2*(f1')^2      [f1-only: ANALYTICALLY computable]
  T2 = f2*(f1')^2          [involves f2]
  T3 = 2*f1*f1'*f2'        [involves f2']
  T4 = -2*f1'*f3'          [involves f3']
  T5 = -(f2')^2            [involves f2']

For T1: compute with mpmath quadosc to full precision (30+ digits).
For T2-T5: compute using scipy ODE + quadosc integration.

Similarly for B5:
S5 = -2f1'f4' - 2f2'f3' + 2f1f1'f3' + f1(f2')^2 + 2f2f1'f2' + f3(f1')^2
     - 2f1^2f1'f2' - 2f1f2(f1')^2 + f1^3(f1')^2

The f1-only terms in S5 are: +f1^3(f1')^2 [ANALYTICALLY computable]
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 40
pi = mp.pi

print("="*78)
print("  B4 DECOMPOSITION: analytical f1-only terms")
print("  dps =", mp.mp.dps)
print("="*78)

def f1(s):
    if s < mp.mpf('1e-30'): return mp.mpf(1)
    return mp.sin(s)/s

def f1p(s):
    if s < mp.mpf('1e-30'): return mp.mpf(0)
    return mp.cos(s)/s - mp.sin(s)/s**2

# ---- T1: -f1^2 * (f1')^2 ----
# This is the "cubic" self-interaction term
# B4_T1 = integral cos(t)*t*(-f1^2*(f1')^2) dt
# A4_T1 = -integral sin(t)*t*(-f1^2*(f1')^2) dt

print("\n  --- T1: -f1^2*(f1')^2 ---")

def integrand_B4_T1(t):
    return mp.cos(t) * t * (-f1(t)**2 * f1p(t)**2)

def integrand_A4_T1(t):
    return -mp.sin(t) * t * (-f1(t)**2 * f1p(t)**2)

B4_T1 = mp.quadosc(integrand_B4_T1, [0, mp.inf], omega=1)
A4_T1 = mp.quadosc(integrand_A4_T1, [0, mp.inf], omega=1)

print(f"  B4_T1 = {mp.nstr(B4_T1, 30)}")
print(f"  A4_T1 = {mp.nstr(A4_T1, 30)}")

# ---- T1 for S5: +f1^3*(f1')^2 ----
print("\n  --- S5 f1-only: f1^3*(f1')^2 ---")

def integrand_B5_f1only(t):
    return mp.cos(t) * t * f1(t)**3 * f1p(t)**2

def integrand_A5_f1only(t):
    return -mp.sin(t) * t * f1(t)**3 * f1p(t)**2

B5_f1only = mp.quadosc(integrand_B5_f1only, [0, mp.inf], omega=1)
A5_f1only = mp.quadosc(integrand_A5_f1only, [0, mp.inf], omega=1)

print(f"  B5_f1only = {mp.nstr(B5_f1only, 30)}")
print(f"  A5_f1only = {mp.nstr(A5_f1only, 30)}")

# ---- Known values for comparison ----
print(f"\n{'='*78}")
print("  KNOWN VALUES for comparison")
print("="*78)

# From r6_c33: B4 ≈ 0.08807035, A4 ≈ -0.10039117 (6 digits from R=10000)
B4_num = mp.mpf('0.088070338')
A4_num = mp.mpf('-0.100391139')
B5_num = mp.mpf('0.006750452')

print(f"  B4 (numerical)   = {mp.nstr(B4_num, 12)}")
print(f"  B4_T1 (f1-only)  = {mp.nstr(B4_T1, 12)}")
print(f"  B4 - B4_T1       = {mp.nstr(B4_num - B4_T1, 12)}")
print(f"  (This is the sum of T2+T3+T4+T5 contributions)")

print(f"\n  A4 (numerical)   = {mp.nstr(A4_num, 12)}")
print(f"  A4_T1 (f1-only)  = {mp.nstr(A4_T1, 12)}")
print(f"  A4 - A4_T1       = {mp.nstr(A4_num - A4_T1, 12)}")

print(f"\n  B5 (numerical)   = {mp.nstr(B5_num, 12)}")
print(f"  B5_f1only        = {mp.nstr(B5_f1only, 12)}")
print(f"  B5 - B5_f1only   = {mp.nstr(B5_num - B5_f1only, 12)}")

# ---- Try PSLQ on T1 terms ----
print(f"\n{'='*78}")
print("  PSLQ on f1-only integrals")
print("="*78)

ln2 = mp.log(2)
ln3 = mp.log(3)

basis_names = ['pi^2', 'pi^4', 'ln2', 'ln3', 'ln2^2', 'ln3^2', 'ln2*ln3',
               'zeta3', 'G', 'pi*ln2', 'pi*ln3', '1']
basis_vals = [pi**2, pi**4, ln2, ln3, ln2**2, ln3**2, ln2*ln3,
              mp.zeta(3), mp.catalan, pi*ln2, pi*ln3, mp.mpf(1)]

for tname, tval in [('B4_T1', B4_T1), ('A4_T1', A4_T1),
                     ('B5_f1only', B5_f1only), ('A5_f1only', A5_f1only)]:
    print(f"\n  --- PSLQ for {tname} = {mp.nstr(tval, 20)} ---")
    basis = [tval] + basis_vals
    names = [tname] + basis_names
    try:
        rel = mp.pslq(basis, maxcoeff=10**8, tol=mp.mpf(10)**(-25))
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

# ---- Also try to evaluate analytically ----
print(f"\n{'='*78}")
print("  ANALYTICAL EVALUATION of T1 integrals")
print("="*78)

# T1: -f1^2*(f1')^2 = -sin^2(t)*(t*cos(t)-sin(t))^2/t^6
# t*T1 = -sin^2(t)*(t*cos(t)-sin(t))^2/t^5
# Expand: (t cos t - sin t)^2 = t^2 cos^2 t - 2t sin t cos t + sin^2 t
# sin^2*(t^2 cos^2-2t sin cos+sin^2) = sin^2 cos^2 t^2 - 2 sin^3 cos t + sin^4
# /t^5: sin^2 cos^2/t^3 - 2 sin^3 cos/t^4 + sin^4/t^5
#
# So integral cos(t) * [-sin^2 cos^2/t^3 + 2 sin^3 cos/t^4 - sin^4/t^5] dt
# = -I1 + 2*I2 - I3 where:
# I1 = integral cos(t)*sin^2(t)*cos^2(t)/t^3 dt
# I2 = integral cos(t)*sin^3(t)*cos(t)/t^4 dt
# I3 = integral cos(t)*sin^4(t)/t^5 dt

# Verify decomposition numerically
I1 = mp.quadosc(lambda t: mp.cos(t)*mp.sin(t)**2*mp.cos(t)**2/t**3, [0, mp.inf], omega=1)
I2 = mp.quadosc(lambda t: mp.cos(t)*mp.sin(t)**3*mp.cos(t)/t**4, [0, mp.inf], omega=1)
I3 = mp.quadosc(lambda t: mp.cos(t)*mp.sin(t)**4/t**5, [0, mp.inf], omega=1)

B4_T1_check = -I1 + 2*I2 - I3
print(f"\n  I1 = {mp.nstr(I1, 30)}")
print(f"  I2 = {mp.nstr(I2, 30)}")
print(f"  I3 = {mp.nstr(I3, 30)}")
print(f"  B4_T1 = -I1 + 2*I2 - I3 = {mp.nstr(B4_T1_check, 30)}")
print(f"  Direct:                    {mp.nstr(B4_T1, 30)}")
print(f"  diff = {mp.nstr(B4_T1_check - B4_T1, 6)}")

# Use trig identities for simplification:
# sin^2 cos^2 = sin^2(2t)/4
# sin^2 cos^3 = sin^2(1+cos(2t))/2... no, let's use product-to-sum
# cos*sin^2*cos^2 = cos * sin^2(2t)/4 = [cos(t)*sin^2(2t)]/4
# sin^2(2t) = (1-cos(4t))/2
# cos(t)*(1-cos(4t))/2 = [cos(t) - cos(t)cos(4t)]/2 = [cos(t) - (cos(3t)+cos(5t))/2]/2
# = [cos(t) - cos(3t)/2 - cos(5t)/2]/2
# = cos(t)/2 - cos(3t)/4 - cos(5t)/4

# I1 = integral [cos(t)/2 - cos(3t)/4 - cos(5t)/4] / t^3 dt
# Can use int cos(at)/t^3 = ?? (divergent at 0!)
# Need to regularize. Actually int cos(at)sin^2(bt)/t^3 is convergent because
# sin^2(bt) ~ b^2 t^2 near 0, making the integrand ~ b^2/t near 0... still divergent!

# Wait, but our original integral IS convergent. Let me check near t=0:
# cos(t)*sin^2(t)*cos^2(t)/t^3 ~ 1*t^2*1/t^3 = 1/t → divergent at 0!
# That can't be right. Let me recheck.

# Actually, t*(-f1^2*(f1')^2) = -sin^2(t)*(t cos t - sin t)^2/t^5
# Near t=0: sin ~ t-t^3/6, t cos t - sin t ~ t(1-t^2/2) - (t-t^3/6) = -t^3/3
# So (t cos t - sin t)^2 ~ t^6/9
# sin^2 ~ t^2
# /t^5: t^2 * t^6/9 / t^5 = t^3/9 → 0. OK, convergent.
# But individually I1 etc may diverge. The cancellation matters.

# Let me try PSLQ on the combined B4_T1 directly
print(f"\n  Note: individual I1, I2, I3 may diverge logarithmically at 0,")
print(f"  but the combination -I1 + 2*I2 - I3 converges (cancellation).")
print(f"  PSLQ result on B4_T1 above is the relevant test.")

# Repeat PSLQ with extended basis including higher-order constants
print(f"\n  Extended PSLQ (with ln^3_2, ln^3_3, pi^2*ln2, ...):")
ext_names = ['pi^2', 'pi^4', 'pi^6', 'ln2', 'ln3', 'ln2^2', 'ln3^2', 'ln2*ln3',
             'zeta3', 'G', 'pi*ln2', 'pi*ln3', 'pi^2*ln2', 'pi^2*ln3', '1']
ext_vals = [pi**2, pi**4, pi**6, ln2, ln3, ln2**2, ln3**2, ln2*ln3,
            mp.zeta(3), mp.catalan, pi*ln2, pi*ln3, pi**2*ln2, pi**2*ln3, mp.mpf(1)]

for tname, tval in [('B4_T1', B4_T1), ('A4_T1', A4_T1)]:
    basis = [tval] + ext_vals
    names = [tname] + ext_names
    try:
        rel = mp.pslq(basis, maxcoeff=10**6, tol=mp.mpf(10)**(-25))
        if rel and rel[0] != 0:
            parts = [f"{c:+d}*{n}" for n,c in zip(names[1:], rel[1:]) if c != 0]
            print(f"    {rel[0]}*{tname} = {', '.join(parts)}")
            check = sum(c*v for c,v in zip(rel, basis))
            print(f"    residual = {mp.nstr(check, 6)}")
        else:
            print(f"    {tname}: no relation")
    except Exception as e:
        print(f"    {tname}: error {e}")

print(f"\n{'='*78}")
