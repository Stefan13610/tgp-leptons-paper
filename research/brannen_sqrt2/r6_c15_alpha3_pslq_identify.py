#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c15_alpha3_pslq_identify.py:

PSLQ integer-relation hunt for alpha_3 against a rich basis of "natural" constants.
Uses the ultra-precise swap-based value alpha_3 = 0.08972222367362532604749467...
(dps=50 from previous session).

Strategy:
1. Direct: try to match alpha_3 = p/q * pi^2 for integer p, q bounded by ~10^5.
2. Extended basis: [alpha_3, pi^2, ln^2(2), ln^2(3), ln(2)ln(3), pi^2·ln(2), pi^2·ln(3), zeta(3), G, chi_2(1/3)]
   and see if PSLQ finds a relation with small integer coefficients.
3. Rational multiplier scan: for each small c, check if c*alpha_3 is a neat rational multiple of pi^2.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 50
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)
zeta3 = mp.zeta(3)
G = mp.catalan

def chi2(x):
    return (mp.polylog(2, x) - mp.polylog(2, -x))/2

chi2_13 = chi2(mp.mpf(1)/3)
chi2_12 = chi2(mp.mpf(1)/2)

# Alpha_3 from the most precise swap-based run to date (r6_c11 at dps=50)
alpha3 = mp.mpf('0.08972222367362532604749467880483973500')

print("="*78)
print("  PSLQ closed-form identification for alpha_3")
print(f"  alpha_3 (numerical, dps=50) = {mp.nstr(alpha3, 40)}")
print("="*78)

# ----- 1. Rational multiple of pi^2 scan -----
print("\n  [1] Simple rational multiple of pi^2:  alpha_3 * N / pi^2 for various N:")
ratio = alpha3 / pi**2
print(f"    alpha_3 / pi^2 = {mp.nstr(ratio, 30)}")
print(f"    1/ratio       = {mp.nstr(1/ratio, 20)}")
# Best near 109.9980...
for N in range(100, 130):
    v = alpha3 * N / pi**2
    if abs(v - round(float(v))) < 0.01:
        dev = v - round(float(v))
        print(f"    N={N}: alpha_3*{N}/pi^2 = {mp.nstr(v, 15)}, dev from nearest int = {mp.nstr(dev, 8)}")

# ----- 2. PSLQ: look for integer rel among a basis -----
print("\n  [2] PSLQ on expanded basis:")
basis_names = ['alpha_3', 'pi^2', 'pi^2*ln2', 'pi^2*ln3', 'ln2^2', 'ln3^2', 'ln2*ln3',
               'zeta(3)', 'Catalan', 'chi2(1/3)', 'chi2(1/2)', '1']
basis = [alpha3, pi**2, pi**2*ln2, pi**2*ln3, ln2**2, ln3**2, ln2*ln3,
         zeta3, G, chi2_13, chi2_12, mp.mpf(1)]

try:
    rel = mp.pslq(basis, maxcoeff=10**8, tol=mp.mpf(10)**(-30))
    if rel:
        print(f"    Found PSLQ relation (max coeff {max(abs(r) for r in rel)}):")
        for name, coeff in zip(basis_names, rel):
            if coeff != 0:
                print(f"       {coeff:>+12d} * {name}")
    else:
        print("    No small-integer relation found.")
except Exception as e:
    print(f"    PSLQ error: {e}")

# ----- 3. Restricted PSLQ: [alpha3, pi^2, ln2, ln3, 1] -----
print("\n  [3] PSLQ restricted to [alpha_3, pi^2, ln2, ln3, 1]:")
b_small = [alpha3, pi**2, ln2, ln3, mp.mpf(1)]
names_small = ['alpha_3', 'pi^2', 'ln2', 'ln3', '1']
try:
    rel = mp.pslq(b_small, maxcoeff=10**10, tol=mp.mpf(10)**(-35))
    if rel:
        for n, c in zip(names_small, rel):
            if c != 0:
                print(f"       {c:>+15d} * {n}")
    else:
        print("    No relation.")
except Exception as e:
    print(f"    PSLQ error: {e}")

# ----- 4. Test: alpha3 vs pi^2/110 - correction -----
print("\n  [4] Correction form: alpha_3 = pi^2/110 + delta")
delta = alpha3 - pi**2/110
print(f"    delta = {mp.nstr(delta, 25)}")
# PSLQ delta
basis_d = [delta, pi**2, ln2, ln3, ln2**2, ln2*ln3, ln3**2, zeta3, G, mp.mpf(1)]
names_d = ['delta', 'pi^2', 'ln2', 'ln3', 'ln2^2', 'ln2*ln3', 'ln3^2', 'zeta3', 'G', '1']
try:
    rel = mp.pslq(basis_d, maxcoeff=10**12, tol=mp.mpf(10)**(-30))
    if rel:
        print(f"    PSLQ for delta:")
        for n, c in zip(names_d, rel):
            if c != 0:
                print(f"       {c:>+20d} * {n}")
    else:
        print("    No relation for delta.")
except Exception as e:
    print(f"    PSLQ error: {e}")

# ----- 5. Alternate rationals near 110 -----
print("\n  [5] Rationals p/q with |p/q - alpha_3/pi^2| < 1e-7, small q:")
target = float(ratio)  # ~0.009090743
# Find close rationals via continued fractions
cf = mp.identify(ratio, ['pi^2'])
print(f"    mpmath identify(ratio, ['pi^2']): {cf}")
# Try identify with more bases
cf2 = mp.identify(alpha3, ['pi**2', 'pi', 'ln(2)', 'ln(3)'])
print(f"    mpmath identify(alpha3, [...]): {cf2}")

print("="*78)
