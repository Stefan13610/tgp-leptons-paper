#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c18_pslq_extended.py:
Extended PSLQ with Clausen Cl_2, inverse tangent Ti_2, Li_3 values.
The integrals involve sin(kt)/t^n * Ci(mt) products which naturally
produce these special values at arguments related to the 1:3 ratio
(from sinc expansion).
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 45
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

# Clausen function Cl_2(theta) = -int_0^theta ln|2sin(t/2)| dt = Im[Li_2(e^{i theta})]
def Cl2(theta):
    return mp.im(mp.polylog(2, mp.expj(theta)))

# Inverse tangent integral Ti_2(x) = int_0^x arctan(t)/t dt = Im[Li_2(ix)]
def Ti2(x):
    return mp.im(mp.polylog(2, 1j*x))

# Higher Clausen Cl_3
def Cl3(theta):
    return mp.re(mp.polylog(3, mp.expj(theta)))

alpha3 = mp.mpf('0.08972222367362532604749467880484')
Pcos = mp.mpf('0.012615939290114711837850')
KcII = mp.mpf('-0.031879037931728580134156')

# Extended basis constants
constants = {
    'pi^2': pi**2,
    'pi^4': pi**4,
    'ln2': ln2,
    'ln3': ln3,
    'ln2^2': ln2**2,
    'ln3^2': ln3**2,
    'ln2*ln3': ln2*ln3,
    'zeta3': mp.zeta(3),
    'Catalan': mp.catalan,
    'Cl2(pi/3)': Cl2(pi/3),
    'Cl2(2pi/3)': Cl2(2*pi/3),
    'Cl2(pi/6)': Cl2(pi/6),
    'Ti2(1/3)': Ti2(mp.mpf(1)/3),
    'Ti2(1/sqrt3)': Ti2(1/mp.sqrt(3)),
    'Li3(1/3)': mp.polylog(3, mp.mpf(1)/3),
    'Li3(-1/3)': mp.polylog(3, mp.mpf(-1)/3),
    'Li3(1/2)': mp.polylog(3, mp.mpf(1)/2),
    'pi*ln2': pi*ln2,
    'pi*ln3': pi*ln3,
    'gamma': mp.euler,
    '1': mp.mpf(1),
}

print("="*78)
print("  EXTENDED PSLQ with Clausen/Ti/Li_3 constants")
print("="*78)
print(f"\n  Basis ({len(constants)} elements):")
for name, val in constants.items():
    print(f"    {name:20s} = {mp.nstr(val, 20)}")

targets = {'alpha_3': alpha3, 'P_cos': Pcos, 'K_c^(II)': KcII}

for tname, tval in targets.items():
    print(f"\n  --- PSLQ for {tname} = {mp.nstr(tval, 25)} ---")
    basis = [tval] + list(constants.values())
    names = [tname] + list(constants.keys())

    for mc in [10**6, 10**8]:
        try:
            rel = mp.pslq(basis, maxcoeff=mc, tol=mp.mpf(10)**(-25))
            if rel and rel[0] != 0:
                parts = []
                for n, c in zip(names[1:], rel[1:]):
                    if c != 0:
                        parts.append(f"{c:+d}*{n}")
                print(f"    maxcoeff={mc}: {rel[0]}*{tname} = {', '.join(parts)}")
                # Verify
                check = sum(c*v for c,v in zip(rel, basis))
                print(f"    residual = {mp.nstr(check, 6)}")
                break
            elif rel:
                print(f"    maxcoeff={mc}: trivial (coeff of target = 0)")
            else:
                print(f"    maxcoeff={mc}: no relation")
        except Exception as e:
            print(f"    maxcoeff={mc}: error: {e}")

print("="*78)
