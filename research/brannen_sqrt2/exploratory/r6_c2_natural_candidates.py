#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c2_natural_candidates.py
============================
Analiza kilku 'naturalnych' kandydatow closed-form dla alpha_3.

Wczesniejszy skan wykryl dopasowanie przy 1.1e-5:
  alpha_3 = ln(3)/48 + ln(3)^2/18
         = (ln3 * (3 + 8 ln3))/144
         = (1-c1)/12 + 8 (1-c1)^2 / 9      (bo ln3 = 4(1-c1))

Ale ODE bias wynosi 2.3e-4, wiec to dopasowanie jest WEWNATRZ szumu.
Kilka innych kandydatow tez miesci sie w tym przedziale.

Cel: znalezc kandydata ktory (a) pasuje lepiej niz szum, i (b) ma
'naturalne' wspolczynniki.

Author: Claudian
Date: 2026-04-16
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import math

ln3 = math.log(3)
ln2 = math.log(2)
c1 = 1 - ln3/4

print("=" * 78)
print("  NATURALNE KANDYDATY dla alpha_3")
print("=" * 78)

alpha3_rich = 0.089939487  # Richardson d -> 0
alpha3_fits = 0.089710     # polynomial fits (all/most points)
alpha3_mid  = (alpha3_rich + alpha3_fits) / 2
print(f"\n  Pomiar numeryczny:")
print(f"    Richardson (d->0): {alpha3_rich:.10f}")
print(f"    poly d^2 fits:     {alpha3_fits:.10f}")
print(f"    srednia:           {alpha3_mid:.10f}")
print(f"    systematyczny bias ~ 2.3e-4")

cands = {
    '(1-c1)/12 + 8(1-c1)^2/9':  (1-c1)/12 + 8*(1-c1)**2/9,
    'ln(3)/48 + ln(3)^2/18':    ln3/48 + ln3**2/18,

    '1/(4*pi)':                 1/(4*math.pi),
    '(pi-1)/24':                (math.pi-1)/24,
    'pi^2/110':                 math.pi**2 / 110,
    'pi^3/346':                 math.pi**3 / 346,
    'pi/35':                    math.pi / 35,

    'c1/8':                     c1/8,
    'c1^2/6':                   c1**2/6,
    'c1^3/4':                   c1**3/4,
    'c1*(1-c1)/2':              c1*(1-c1)/2,
    '(1-c1)/3':                 (1-c1)/3,
    '(1-c1)^2*5/4':             (1-c1)**2*5/4,
    '1/2 - 3*c1/5':             0.5 - 3*c1/5,

    'ln(2)/8':                  ln2/8,
    'ln(2)^2/5':                ln2**2/5,
    'ln(3)/12':                 ln3/12,
    'ln(3)^2/16':               ln3**2/16,

    '1/12 + 1/120':             1/12 + 1/120,
    '1/12 + 1/180':             1/12 + 1/180,
    '1/12 + c1^2/100':          1/12 + c1**2/100,

    'c1^2 / 8':                 c1**2/8,
    '(c1/2)^2':                 (c1/2)**2,

    '1/3 - c1/4':               1/3 - c1/4,
    '(2 - ln3)/10':             (2 - ln3)/10,
    '(5 - 3*ln3)/18':           (5 - 3*ln3)/18,
    '(3 - 2*ln3)/(4*pi)':       (3 - 2*ln3)/(4*math.pi),
}

def eval_cand(name, val):
    d_rich = val - alpha3_rich
    d_fit = val - alpha3_fits
    d_mid = val - alpha3_mid
    best = min(abs(d_rich), abs(d_fit), abs(d_mid))
    return (name, val, d_rich, d_fit, best)

results = [eval_cand(n, v) for n, v in cands.items()]
results.sort(key=lambda x: x[4])

print(f"\n  {'Kandydat':30s}  {'Wartosc':>14s}  {'d(rich)':>11s}  {'d(fit)':>11s}  {'best':>10s}")
print(f"  {'-'*30}  {'-'*14}  {'-'*11}  {'-'*11}  {'-'*10}")
for name, val, d_rich, d_fit, best in results[:20]:
    inside_noise = " <--" if best < 2.3e-4 else ""
    print(f"  {name:30s}  {val:14.10f}  {d_rich:+11.2e}  {d_fit:+11.2e}  {best:10.2e}{inside_noise}")

inside = [r for r in results if r[4] < 2.3e-4]
print(f"\n  KANDYDATOW wewnatrz bias 2.3e-4: {len(inside)}")
print(f"  Zadnego nie mozna wyrozniac bez poprawy precyzji.")
print(f"\n  Aby zidentyfikowac closed-form, potrzeba redukcji bias o rzad ~10x.")
