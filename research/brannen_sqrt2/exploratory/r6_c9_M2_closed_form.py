#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c9_M2_closed_form.py: Analytical derivation of M2(a,b) = int sin(a*t) Ci(b*t)/t dt.

Derivation (differentiating under integral):
    F(a) := int_0^inf sin(a t) Ci(b t)/t dt    (b > 0 fixed)
    F'(a) = int cos(a t) Ci(b t) dt
          = -(1/(2a))[sgn(a+b)·pi/2 + sgn(a-b)·pi/2]    [IBP + Dirichlet]

For a, b > 0:
    - If a < b:  F'(a) = 0           =>  F(a) = F(0) = 0
    - If a > b:  F'(a) = -pi/(2a)    =>  F(a) = -(pi/2)ln(a/b)

So M2(a,b) = -(pi/2) ln(a/b) * Heaviside(a - b).
  Symbolically: M2(a,b) = -(pi/2) ln(max(a,b)/b) (with indicator).

Similarly for M3(a,b) = int cos(a*t) Si(b*t)/t dt and M4(a,b) = int cos(a*t) Ci(b*t)/t dt:
These are singular at t=0 for some combinations; need regularization.

Tests: M2(1,1) = 0, M2(2,1) = -(pi/2)ln(2), M2(3,1) = -(pi/2)ln(3), M2(1,3) = 0.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp
mp.mp.dps = 30

pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)
Ci = mp.ci
Si = mp.si

print("="*72)
print("  M2(a,b) = int sin(a t) Ci(b t)/t dt    CLOSED FORM")
print("    M2(a,b) = -(pi/2) ln(a/b) if a > b, else 0.")
print("="*72)

def M2_th(a, b):
    if a > b:
        return -pi/2 * mp.log(mp.mpf(a)/b)
    else:
        return mp.mpf(0)

def M2_num(a, b):
    return mp.quadosc(lambda t: mp.sin(a*t) * Ci(b*t)/t, [mp.mpf('1e-30'), mp.inf], period=2*pi/max(a,b))

test_cases = [(1,1), (1,2), (2,1), (1,3), (3,1), (2,3), (3,2), (4,1), (1,4)]
print("\n  Tests:")
for a, b in test_cases:
    th = M2_th(a, b)
    num = M2_num(a, b)
    print(f"    M2({a},{b}): th = {mp.nstr(th, 20):>22}  num = {mp.nstr(num, 20):>22}  diff = {mp.nstr(th-num, 4)}")

# ==== M3(a,b) = int cos(a t) Si(b t)/t dt ====
# Derivation: F(b) := int cos(a t) Si(b t)/t dt (a > 0 fixed)
# F'(b) = int cos(a t) * [d/db Si(b t)]/t dt = int cos(a t) sin(b t)/(b t) * b dt = int cos(a t) sin(b t)/t dt
# Wait: d/db Si(b t) = sin(b t) * t / (b t) * ?? No: d/dx Si(x) = sin(x)/x, so d/db Si(b t) = sin(bt)/bt * t = sin(bt)/b
# So F'(b) = int cos(a t) sin(b t)/(b t) dt.
# Hmm that's int cos(at) sin(bt)/(bt) dt = (1/b) int cos(at) sin(bt)/t dt.
# And int cos(a t) sin(b t)/t dt = (1/2)[Si-style]:
#   sin(bt)cos(at) = [sin((b+a)t) + sin((b-a)t)]/2
#   int sin((b+a)t)/t dt = pi/2·sgn(b+a) = pi/2  (for a,b > 0)
#   int sin((b-a)t)/t dt = pi/2·sgn(b-a)
# So int cos(at)sin(bt)/t dt = pi/4·[1 + sgn(b-a)]
#   = pi/2 for b > a,  pi/4 for b = a,  0 for b < a.
# Therefore F'(b) = (1/b) * {pi/2 for b>a, pi/4 at b=a, 0 for b<a}.
#   = pi/(2b) for b > a,  0 for b < a.
# Integrating: F(b) = F(a) + int_a^b pi/(2b') db' = F(a) + (pi/2) ln(b/a) for b > a.
# But we also know F(b)|_{b=0} = 0 (since Si(0) = 0).
# For b < a, F'(b) = 0, so F(b) = F(0) = 0. Thus F(a) = lim_{b -> a-} F(b) = 0.
# So F(b) = (pi/2) ln(b/a) for b > a, F(b) = 0 for b < a.
# That is M3(a,b) = -(pi/2) ln(a/b) ... wait no, = (pi/2) ln(b/a) for b > a, else 0.

print(f"\n{'='*72}")
print(f"  M3(a,b) = int cos(a t) Si(b t)/t dt    CLOSED FORM")
print(f"    M3(a,b) = (pi/2) ln(b/a) if b > a, else 0.")
print(f"{'='*72}")

def M3_th(a, b):
    if b > a:
        return pi/2 * mp.log(mp.mpf(b)/a)
    else:
        return mp.mpf(0)

def M3_num(a, b):
    return mp.quadosc(lambda t: mp.cos(a*t) * Si(b*t)/t, [mp.mpf('1e-30'), mp.inf], period=2*pi/max(a,b))

print("\n  Tests:")
for a, b in test_cases:
    th = M3_th(a, b)
    num = M3_num(a, b)
    print(f"    M3({a},{b}): th = {mp.nstr(th, 20):>22}  num = {mp.nstr(num, 20):>22}  diff = {mp.nstr(th-num, 4)}")

print(f"\n{'='*72}")
print(f"  NOW M4(a,b) = int cos(a t) Ci(b t)/t dt -- likely DIVERGENT at t=0 (Ci~ln)")
print(f"{'='*72}")
# cos(at)Ci(bt)/t near t=0: ~1·(gamma + ln(bt))/t → diverges. So M4 is only PV-regularized.
# Consider instead M4-like integrals in differences: M4(a, b_1) - M4(a, b_2) = -(pi/2)·... hmm.
# In our A_i we only ever see Ci(t) - Ci(3t), which is finite at t=0.
# So we can define M4_diff(a; b1, b2) = int cos(a t) [Ci(b1 t) - Ci(b2 t)]/t dt, which is finite.
# d/db Ci(bt) = cos(bt)·t / bt · ... no: d/dx Ci(x) = cos(x)/x, so d/db Ci(bt) = cos(bt)/b.
# So M4(a,b) - M4(a,b_0) = int cos(a t) int_{b_0}^b cos(c t)/c dc / t dt
#                       = int_{b_0}^b (dc/c) int cos(a t) cos(c t)/t dt
# int cos(a t) cos(c t)/t dt: diverges at 0. Need regularization.
# Actually cos(at)cos(ct) = [cos((a+c)t) + cos((a-c)t)]/2.
# ∫_0^∞ cos(kt)/t dt is divergent (Ci divergent).

# Alternative: for difference of Ci,
#   Ci(b1 t) - Ci(b2 t) = -int_{b1 t}^{b2 t} cos(u)/u du = -ln(b2/b1) (leading) at t→0.
# Actually Ci(bt) = gamma + ln(bt) + Integral, so Ci(b1 t) - Ci(b2 t) = ln(b1/b2) + O(t²).
# Converges at 0. Good.

# So let's derive: Let G(b) := int cos(a t) Ci(b t)/t dt (PV or with cutoff).
# dG/db = int cos(a t) cos(b t)/(b t) * t? No: d/db Ci(bt) = cos(bt)/b (treating bt as x, dCi/dx = cos(x)/x, d/db = t·cos(bt)/(bt) = cos(bt)/b).
# So dG/db = (1/b) int cos(a t) cos(b t)/t dt.
# int cos(at)cos(bt)/t dt = (1/2) int [cos((a+b)t) + cos((a-b)t)]/t dt — divergent.
# But the DIFFERENCE G(b1) - G(b2) with b1 ≠ b2 might be finite via regularized evaluation.

# Actually using standard tables (Gradshteyn 4.421):
# int_0^∞ cos(ax) [Ci(bx) - Ci(cx)] dx = pi/(2a) * [Heaviside(a-b) - Heaviside(a-c)] ?
# Or with /t factor: different formulas.

# Simpler: for our purposes, we need int cos(at)·(Ci(3t)-Ci(t))/t dt. Let's derive.
# By IBP: u = Ci(3t) - Ci(t), du = [cos(3t) - cos(t)]/t dt.
# dv = cos(at)/t dt — divergent antiderivative. Instead:
# Switch: d/dt[Ci(3t) - Ci(t)] = cos(3t)/t - cos(t)/t = [cos(3t)-cos(t)]/t.
# IBP: u = [Ci(3t)-Ci(t)]/t? Messy.

# Try Parseval: Ci(bt) - Ci(at) = -int_a^b cos(ct)/c dc. So:
# int cos(pt)(Ci(b t) - Ci(a t))/t dt = -int_a^b (dc/c) int cos(p t) cos(c t)/t dt
# But inner diverges.

# Back to F(b) for M3-type: we had a clean derivation. Let me handle M4 similarly with difference.

print("  For M4(a, b_1) - M4(a, b_2) = int cos(a t) [Ci(b_1 t) - Ci(b_2 t)]/t dt (FINITE):")
print("  Approach: differentiate w.r.t. a (not b):")
print("    G(a) := int cos(a t) [Ci(b_1 t) - Ci(b_2 t)]/t dt")
print("    dG/da = -int sin(a t) [Ci(b_1 t) - Ci(b_2 t)] dt")
print("    Still requires explicit evaluation. Use Gradshteyn 4.421 or similar.")

# Numerical verification:
def M4diff_num(a, b1, b2):
    return mp.quadosc(lambda t: mp.cos(a*t) * (Ci(b1*t) - Ci(b2*t))/t,
                      [mp.mpf('1e-30'), mp.inf], period=2*pi/max(a,b1,b2))

print("\n  Numerical tests of M4diff(a, b1, b2):")
for a, b1, b2 in [(1, 3, 1), (2, 3, 1), (3, 3, 1), (4, 3, 1), (5, 3, 1)]:
    val = M4diff_num(a, b1, b2)
    print(f"    M4diff({a}, {b1}, {b2}) = {mp.nstr(val, 20)}")

# Look for pattern: compare M4diff with something like -(pi/2)ln(...)
# For a=1, b1=3, b2=1: a < b2, so "both outside" — might be -ln(3)?
# a between b2 and b1: transition zone.
# a > b1: ?

# Intuition via F(b) ansatz: M4(a,b_1) - M4(a,b_2) = -(pi/2)[H(a-b_2) ln(a/b_2)- H(a-b_1) ln(a/b_1)] ???

# Let me try a guessed formula:
# M4diff(a, b1, b2) = ??? ln(b1/b2) * piecewise function of a vs b's.

# Check a=1, b1=3, b2=1: value ≈ ?
# Guess: for a < b2 < b1, M4diff = -ln(b1/b2)·(pi/???)... or 0.
# Actually from M3 derivation: M3(a,b) = (pi/2) ln(b/a) for b > a, 0 for b < a.
# By analogy, M4(a,b) "should" be something like (pi/2) ln(max(a,b)/b) ...? But M4 has divergence.

print("="*72)
