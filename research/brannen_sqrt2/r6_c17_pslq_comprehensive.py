#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c17_pslq_comprehensive.py:
Comprehensive PSLQ identification for alpha_3, P_cos, K_c^(II), and individual A_i.
Uses the CONFIRMED 30-digit values from the dps=60 convergence study.

Basis: pi^2, pi^4, ln(2), ln(3), ln^2(2), ln^2(3), ln(2)ln(3), zeta(3),
       Catalan G, chi_2(1/3), chi_2(1/2), pi*ln(2), pi*ln(3), euler_gamma, 1.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 45
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)
gamma = mp.euler
zeta3 = mp.zeta(3)
G = mp.catalan

def chi2(x):
    return (mp.polylog(2, x) - mp.polylog(2, -x))/2

chi2_13 = chi2(mp.mpf(1)/3)
chi2_12 = chi2(mp.mpf(1)/2)

# Confirmed values (30 digits, from dps=60 run)
A1 = mp.mpf('0.0703145356797345181974810338624')
A2 = mp.mpf('0.0291881503799267524695228131641')
A3 = mp.mpf('0.00454851496514307720777933676503')
A4 = mp.mpf('-0.0630751331627896133250687447822')

KcI = (ln2 - 1) / 6
KcII = -A1 - A2 + A3 - A4
Pcos = KcI - 2*KcII
alpha3 = pi**2/128 + Pcos

print("="*78)
print("  COMPREHENSIVE PSLQ IDENTIFICATION (confirmed 30 digits)")
print("="*78)
print(f"  alpha_3   = {mp.nstr(alpha3, 32)}")
print(f"  P_cos     = {mp.nstr(Pcos, 32)}")
print(f"  K_c^(II)  = {mp.nstr(KcII, 32)}")
print(f"  A_1       = {mp.nstr(A1, 32)}")
print(f"  A_2       = {mp.nstr(A2, 32)}")
print(f"  A_3       = {mp.nstr(A3, 32)}")
print(f"  A_4       = {mp.nstr(A4, 32)}")

# ---- Define PSLQ basis ----
basis_full_names = [
    'X', 'pi^2', 'pi^4', 'ln2', 'ln3', 'ln2^2', 'ln3^2', 'ln2*ln3',
    'zeta3', 'G', 'chi2(1/3)', 'chi2(1/2)', 'pi*ln2', 'pi*ln3', 'gamma', '1'
]
basis_full_vals = [
    None,  # placeholder for target
    pi**2, pi**4, ln2, ln3, ln2**2, ln3**2, ln2*ln3,
    zeta3, G, chi2_13, chi2_12, pi*ln2, pi*ln3, gamma, mp.mpf(1)
]

def try_pslq(target, target_name, max_coeff=10**8, tol_exp=-25):
    basis = [target] + basis_full_vals[1:]
    names = [target_name] + basis_full_names[1:]
    try:
        rel = mp.pslq(basis, maxcoeff=max_coeff, tol=mp.mpf(10)**tol_exp)
        if rel:
            coeff_tgt = rel[0]
            if coeff_tgt == 0:
                print(f"    [TRIVIAL: coefficient of {target_name} is 0]")
                return None
            parts = []
            for n, c in zip(names[1:], rel[1:]):
                if c != 0:
                    parts.append(f"{c:+d}*{n}")
            rhs = ", ".join(parts) if parts else "0"
            print(f"    {coeff_tgt}*{target_name} = {rhs}")
            # Verify
            val = sum(c * v for c, v in zip(rel, basis))
            print(f"    residual = {mp.nstr(val, 8)}")
            return rel
        else:
            print(f"    No relation found for {target_name} (maxcoeff={max_coeff})")
            return None
    except Exception as e:
        print(f"    PSLQ error for {target_name}: {e}")
        return None

# ---- Run PSLQ on each target ----
targets = {
    'alpha_3': alpha3,
    'P_cos': Pcos,
    'K_c^(II)': KcII,
    'A_1': A1,
    'A_2': A2,
    'A_3': A3,
    'A_4': A4,
}

for tname, tval in targets.items():
    print(f"\n  --- PSLQ for {tname} = {mp.nstr(tval, 20)} ---")
    # Try with increasing maxcoeff
    for mc in [10**5, 10**8, 10**10]:
        result = try_pslq(tval, tname, max_coeff=mc)
        if result:
            break

# ---- Restricted bases ----
print(f"\n{'='*78}")
print(f"  RESTRICTED PSLQ (smaller bases)")
print(f"{'='*78}")

# Basis: [X, pi^2, ln2, ln3, 1]
for tname, tval in targets.items():
    basis = [tval, pi**2, ln2, ln3, mp.mpf(1)]
    names = [tname, 'pi^2', 'ln2', 'ln3', '1']
    try:
        rel = mp.pslq(basis, maxcoeff=10**12, tol=mp.mpf(10)**(-30))
        if rel and rel[0] != 0:
            parts = []
            for n, c in zip(names[1:], rel[1:]):
                if c != 0:
                    parts.append(f"{c:+d}*{n}")
            print(f"  {rel[0]}*{tname} = {', '.join(parts)}")
        else:
            print(f"  {tname}: no simple relation with [pi^2, ln2, ln3, 1]")
    except:
        print(f"  {tname}: PSLQ failed")

# ---- Special: test K_c^(II) against (ln2-1)/12 + c*pi^2 ----
print(f"\n{'='*78}")
print(f"  K_c^(II) near (ln2-1)/12 + c*pi^2:")
print(f"{'='*78}")
correction = KcII - (ln2-1)/12
print(f"  K_c^(II) - (ln2-1)/12 = {mp.nstr(correction, 25)}")
ratio = correction / pi**2
print(f"  (K_c^(II) - (ln2-1)/12) / pi^2 = {mp.nstr(ratio, 25)}")
# Check if ratio is a simple fraction
print(f"  1/ratio = {mp.nstr(1/ratio, 20)}")
# Closest simple fractions
for d in range(1, 30):
    for n in range(-30, 0):
        if abs(float(ratio) - n/d) < 1e-5:
            print(f"    candidate: {n}/{d} = {n/d:.10f}, diff = {float(ratio) - n/d:.2e}")

# ---- Try: alpha_3 = rational? ----
print(f"\n  alpha_3 as rational (continued fraction):")
cf = mp.identify(alpha3)
print(f"    mpmath.identify: {cf}")

print(f"\n  P_cos as rational (continued fraction):")
cf2 = mp.identify(Pcos)
print(f"    mpmath.identify: {cf2}")

# ---- Try specific combinations ----
print(f"\n{'='*78}")
print(f"  SPECIFIC CLOSED-FORM TESTS:")
print(f"{'='*78}")

tests = {
    'pi^2/110': pi**2/110,
    '(ln2-1)/6 + pi^2/64 - (ln3)^2/8': (ln2-1)/6 + pi**2/64 - ln3**2/8,
    '(pi^2-1)/110': (pi**2-1)/110,
    'pi^2/128 + (ln2-1)/6 + ln3^2/16': pi**2/128 + (ln2-1)/6 + ln3**2/16,
    '3*ln3^2/32 - ln2/4 + pi^2/128 + 1/6': 3*ln3**2/32 - ln2/4 + pi**2/128 + mp.mpf(1)/6,
    'pi^2/128 + (3*pi^2 - 32)/7040': pi**2/128 + (3*pi**2 - 32)/7040,
}
for name, val in tests.items():
    diff = alpha3 - val
    print(f"  {name:45s} = {mp.nstr(val, 20)}  diff = {mp.nstr(diff, 8)}")

print("="*78)
