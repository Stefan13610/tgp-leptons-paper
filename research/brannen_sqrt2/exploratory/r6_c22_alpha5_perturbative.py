#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c22_alpha5_perturbative.py:
High-precision computation of alpha_5 via perturbation expansion to O(delta^4).

Expansion: g = 1 + delta*f1 + delta^2*f2 + delta^3*f3 + delta^4*f4 + ...

Each fn satisfies: fn'' + (2/r)fn' + fn = Source_n(f1,...,f_{n-1})

The source terms come from expanding (1/g)(g')^2 to each order:

O(delta^1): f1'' + (2/r)f1' + f1 = 0
O(delta^2): f2'' + (2/r)f2' + f2 = -(f1')^2
O(delta^3): f3'' + (2/r)f3' + f3 = -(2f1'f2' - f1(f1')^2)
O(delta^4): f4'' + (2/r)f4' + f4 = -((f2')^2 + 2f1'f3' - 2f1f1'f2' - f2(f1')^2 + f1^2(f1')^2)

Tails: fn(r) ~ (A_n*cos(r) + B_n*sin(r))/r as r->inf
So r*fn ~ A_n*cos(r) + B_n*sin(r)

The full tail: r*(g-1) = delta*(r*f1) + delta^2*(r*f2) + ...

S_total = delta*B1 + delta^2*B2 + delta^3*B3 + delta^4*B4 + ...
C_total = delta*A1 + delta^2*A2 + delta^3*A3 + delta^4*A4 + ...

(using A,B for cos,sin amplitudes of r*fn)

A_tail = sqrt(S^2 + C^2)
eta = A_tail/|delta| = sqrt((B1+dB2+d^2*B3+d^3*B4+...)^2 + (A1+dA2+d^2*A3+d^3*A4+...)^2)

For f1 = sin(r)/r: r*f1 = sin(r), so B1=1, A1=0.

eta_sym (symmetrized):
eta_sym = 1 + alpha_3*d^2 + alpha_5*d^4 + ...

From the expansion:
alpha_3 = B2 + A2^2/2  (= I_cos + I_sin^2/2 = P_cos + pi^2/128... wait)

Actually let me be more careful. We use un = r*fn, then:
u1 = sin(r), so B1=1, A1=0.

From the h-equation (f2), we know:
u2 -> I_cos*sin(r) - I_sin*cos(r) as r->inf
So B2 = I_cos = 1/2 - ln3/8, A2 = -I_sin = pi/8.

From the p-equation (f3), we know:
u3 -> P_cos*sin(r) - P_sin*cos(r) as r->inf
So B3 = P_cos, A3 = -P_sin.

From the q-equation (f4):
u4 -> Q_cos*sin(r) - Q_sin*cos(r) as r->inf
So B4 = Q_cos, A4 = -Q_sin.

Now:
S = delta*(1 + delta*B2 + delta^2*B3 + delta^3*B4 + ...)
C = delta*(delta*A2 + delta^2*A3 + delta^3*A4 + ...)

eta = sqrt(S^2+C^2)/|delta|

For eta_sym = (eta(+d) + eta(-d))/2:
Note S(+d) = d(1+dB2+d^2*B3+d^3*B4), S(-d) = -d(1-dB2+d^2*B3-d^3*B4)
     C(+d) = d(dA2+d^2*A3+d^3*A4),    C(-d) = -d(-dA2+d^2*A3-d^3*A4)

A_tail(+d) = d*sqrt((1+dB2+d^2*B3+d^3*B4)^2 + (dA2+d^2*A3+d^3*A4)^2)
A_tail(-d) = d*sqrt((1-dB2+d^2*B3-d^3*B4)^2 + (-dA2+d^2*A3-d^3*A4)^2)
            = d*sqrt((1-dB2+d^2*B3-d^3*B4)^2 + (dA2-d^2*A3+d^3*A4)^2)

eta(+d) = sqrt((1+x)^2 + y^2) where x = dB2+d^2*B3+d^3*B4, y = dA2+d^2*A3+d^3*A4
eta(-d) = sqrt((1-x'+...)^2 + y'^2) similarly

Let me just compute this symbolically to 4th order in d:

Let S = 1 + a*d + b*d^2 + c*d^3 (where a=B2, b=B3, c=B4)
    C = p*d + q*d^2 + r*d^3   (where p=A2, q=A3, r=A4)

eta = sqrt(S^2 + C^2) = sqrt(1 + 2ad + (a^2+2b+p^2)d^2 + (2ab+2c+2pq)d^3
     + (b^2+2ac+q^2+2pr)d^4 + ...)

= 1 + ad + [(a^2+2b+p^2)/2 - a^2/2]d^2 + ...
= 1 + ad + [b + p^2/2]d^2 + [ab + c + pq - a(b+p^2/2)]d^3
  + ...

For eta_sym, only even powers survive:
eta_sym = 1 + [b + p^2/2]d^2 + [???]d^4

alpha_3 = b + p^2/2 = B3 + A2^2/2

Wait, that gives alpha_3 = B3 + A2^2/2 = P_cos + (pi/8)^2/2 = P_cos + pi^2/128.
YES, that matches!

For alpha_5, I need the d^4 coefficient of eta_sym. Let me expand more carefully:

eta^2 = S^2 + C^2
     = (1+ad+bd^2+cd^3)^2 + (pd+qd^2+rd^3)^2
     = 1 + 2ad + (a^2+2b)d^2 + (2ab+2c)d^3 + (b^2+2ac)d^4
       + p^2*d^2 + 2pq*d^3 + (q^2+2pr)*d^4
     = 1 + 2ad + (a^2+2b+p^2)d^2 + (2ab+2c+2pq)d^3
       + (b^2+2ac+q^2+2pr)d^4 + O(d^5)

Let E2 = a^2+2b+p^2, E3 = 2ab+2c+2pq, E4 = b^2+2ac+q^2+2pr

eta = sqrt(1 + 2ad + E2*d^2 + E3*d^3 + E4*d^4)
    = (1 + X)^{1/2} where X = 2ad + E2*d^2 + E3*d^3 + E4*d^4

(1+X)^{1/2} = 1 + X/2 - X^2/8 + X^3/16 - 5X^4/128 + ...

X = 2ad + E2*d^2 + E3*d^3 + E4*d^4
X^2 = 4a^2*d^2 + 4aE2*d^3 + (E2^2+4aE3)*d^4 + ...
X^3 = 8a^3*d^3 + 12a^2*E2*d^4 + ...
X^4 = 16a^4*d^4 + ...

So eta = 1 + [ad + E2*d^2/2 + E3*d^3/2 + E4*d^4/2]
           - [a^2*d^2/2 + aE2*d^3/2 + (E2^2+4aE3)*d^4/8]
           + [a^3*d^3/2 + 3a^2*E2*d^4/4]
           - [5a^4*d^4/8]
         + O(d^5)

Collecting by power of d:
d^1: a
d^2: E2/2 - a^2/2 = (a^2+2b+p^2)/2 - a^2/2 = b + p^2/2   ✓ (= alpha_3)
d^3: E3/2 - aE2/2 + a^3/2
   = (ab+c+pq) - a(a^2+2b+p^2)/2 + a^3/2
   = ab+c+pq - a^3/2 - ab - ap^2/2 + a^3/2
   = c + pq - ap^2/2
d^4: E4/2 - (E2^2+4aE3)/8 + 3a^2*E2/4 - 5a^4/8
   Let me compute each:
   E4/2 = (b^2+2ac+q^2+2pr)/2
   E2^2 = (a^2+2b+p^2)^2 = a^4+4a^2*b+4b^2+2a^2*p^2+4bp^2+p^4
   4aE3 = 4a(2ab+2c+2pq) = 8a^2*b+8ac+8apq
   (E2^2+4aE3)/8 = (a^4+4a^2*b+4b^2+2a^2*p^2+4bp^2+p^4+8a^2*b+8ac+8apq)/8
                  = (a^4+12a^2*b+4b^2+2a^2*p^2+4bp^2+p^4+8ac+8apq)/8
   3a^2*E2/4 = 3a^2*(a^2+2b+p^2)/4 = (3a^4+6a^2*b+3a^2*p^2)/4
   5a^4/8

   alpha_5 = (b^2+2ac+q^2+2pr)/2
           - (a^4+12a^2*b+4b^2+2a^2*p^2+4bp^2+p^4+8ac+8apq)/8
           + (3a^4+6a^2*b+3a^2*p^2)/4
           - 5a^4/8

   = b^2/2 + ac + q^2/2 + pr
     - a^4/8 - 3a^2*b/2 - b^2/2 - a^2*p^2/4 - bp^2/2 - p^4/8 - ac - apq
     + 3a^4/4 + 3a^2*b/2 + 3a^2*p^2/4
     - 5a^4/8

   = q^2/2 + pr - apq + a^2*p^2/2 - bp^2/2 - p^4/8

   Simplify:
   alpha_5 = q^2/2 + pr - apq + p^2*(a^2-b)/2 - p^4/8

   With a=B2, b=B3, c=B4, p=A2, q=A3, r=A4:
   alpha_5 = A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8

   This can be verified! We know:
   B2 = I_cos = 1/2 - ln3/8
   A2 = -I_sin = pi/8
   B3 = P_cos (known to 30 digits)
   A3 = -P_sin (need to compute)
   B4 = Q_cos (need to compute)
   A4 = -Q_sin (need to compute)

   We need P_sin, Q_cos, Q_sin from perturbation theory.

Strategy:
1. Compute f3 (= p) numerically at high precision, extract P_cos and P_sin
   (P_sin should be extractable from the same code that gave P_cos)
2. Set up and solve f4 equation, extract Q_cos and Q_sin
3. Combine via the formula above to get alpha_5
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 40
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

print("="*78)
print("  ALPHA_5 via perturbation theory O(delta^4)")
print("  Working precision: dps =", mp.mp.dps)
print("="*78)

# Known tail amplitudes (from previous computations)
B2 = mp.mpf(1)/2 - ln3/8                # I_cos (proven)
A2 = pi/8                                # -I_sin = pi/8 (proven)
B3_known = mp.mpf('0.012615939290114711837850')  # P_cos
# A3 = -P_sin: needs computation

print(f"\n  Known quantities:")
print(f"    B2 = I_cos = 1/2 - ln3/8 = {mp.nstr(B2, 25)}")
print(f"    A2 = pi/8 = {mp.nstr(A2, 25)}")
print(f"    B3 = P_cos = {mp.nstr(B3_known, 25)}")

# --- Step 1: Compute P_sin from the p-equation source ---
# p'' + (2/r)p' + p = -(2f1'f2' - f1(f1')^2) =: -R(r)
# v = rp: v'' + v = -rR(r)
# v -> P_cos*sin(r) - P_sin*cos(r)
#
# P_cos = int_0^inf cos(s) * s * R(s) ds   (already known)
# P_sin = -int_0^inf sin(s) * s * R(s) ds
# (sign convention from variation of parameters with Green's function)
#
# Actually from variation of parameters:
# v(r) = sin(r)*int_0^r cos(s)*S(s)ds - cos(r)*int_0^r sin(s)*S(s)ds
# where S(s) = -s*R(s)
# So as r->inf:
# v -> sin(r) * [-int_0^inf cos(s)*s*R(s)ds] - cos(r) * [-int_0^inf sin(s)*s*R(s)ds]
# Wait, sign: S(s) = rhs of v''+v = rhs, where rhs = -s*R(s)
# Actually the ODE is v'' + v = -rR(r), so the source is -rR(r).
# v(r) = sin(r)*int_0^r cos(s)*(-s*R(s))ds - cos(r)*int_0^r sin(s)*(-s*R(s))ds
# v -> -sin(r)*int_0^inf cos(s)*s*R(s)ds + cos(r)*int_0^inf sin(s)*s*R(s)ds
# So: P_cos_coeff_of_sin = -int cos*sR ds
# P_sin_coeff_of_cos = +int sin*sR ds (but with minus sign in convention)
#
# Let me use the Green's function form directly:
# v(r) = int_0^r sin(r-s) * Source(s) ds = int_0^r sin(r-s)*(-s*R(s))ds
# = sin(r)*int_0^r cos(s)*(-sR(s))ds - cos(r)*int_0^r sin(s)*(-sR(s))ds
#
# Hmm, I'm getting confused with signs. Let me compute both integrals and
# check P_cos against the known value to fix the sign convention.

# f1(r) = sin(r)/r
# f1'(r) = cos(r)/r - sin(r)/r^2

# For h (=f2), we solved u=r*h:
# u'' + u = -r*(f1')^2
# u -> sin(r)*J_c - cos(r)*J_s  (from the convention in the previous code)
# where J_c = int_0^inf cos(s)*[-s*(f1')^2]ds = 1/2 - ln3/8 = B2  ✓
# and   J_s = int_0^inf sin(s)*[-s*(f1')^2]ds = pi/8  ✓
# Wait, but J_s was defined as int sin(s)*(-s*(f1')^2)ds = pi/8
# And A2 = pi/8 is the coefficient of -cos(r), matching u -> ...J_c*sin(r) - J_s*cos(r)
# But then the cos coefficient is -J_s = -pi/8, so A2 = -(-J_s)... hmm.
# Let me re-check: u2 = r*f2, and u2 ~ I_cos*sin(r) - I_sin*cos(r)
# I_cos = B2 = J_c = 1/2-ln3/8
# I_sin = -A2 = -pi/8? No...
# Actually I_sin = -pi/8 was from the original computation.
# And A2 = -I_sin = pi/8.
# So u2 ~ I_cos*sin(r) - I_sin*cos(r) = B2*sin(r) + A2*cos(r)
# Wait, I_sin = -pi/8, so -I_sin = pi/8 = A2.
# u2 ~ I_cos*sin(r) - I_sin*cos(r) = I_cos*sin + (pi/8)*cos
# So u2 cos-amplitude = pi/8 = A2. ✓

# For the p (=f3) equation, the source is:
# R(r) = 2*f1'*f2' - f1*(f1')^2
# And the Green's function gives:
# v = r*f3
# v(r) = int_0^r sin(r-s)*(-s*R(s)) ds
# As r->inf: v -> sin(r)*int_0^inf cos(s)*(-sR(s))ds - cos(r)*int_0^inf sin(s)*(-sR(s))ds
# So: sin-amplitude of v = -int cos(s)*s*R(s) ds = P_cos (using same convention as J_c)
# cos-amplitude of v = int sin(s)*s*R(s) ds = -(P_sin) where P_sin = -int sin*sR ds

# OK, let me just use mpmath quadosc to compute both P_cos and P_sin.

# First need f2' = u'/r - u/r^2 where u = r*f2
# For u, we need to evaluate via integral...

# Actually, let me use a different strategy: compute f2(r) and f2'(r) via
# the Green's function integral at arbitrary r, then use that to evaluate
# the P_cos and P_sin integrals.

# Or better: use the swap approach. f2(r) = u(r)/r where
# u(r) = int_0^r sin(r-s)*(-s*(f1'(s))^2) ds

# For the outer integral, I need:
# int_0^inf cos(s) * s * R(s) ds and int_0^inf sin(s) * s * R(s) ds
# where R(s) = 2*f1'(s)*f2'(s) - f1(s)*(f1'(s))^2

# f1(s) = sin(s)/s
# f1'(s) = cos(s)/s - sin(s)/s^2 = (s*cos(s) - sin(s))/s^2

# f2(s) = u(s)/s where u(s) = int_0^s sin(s-t)*(-t*(f1'(t))^2) dt
# f2'(s) = (u'(s)*s - u(s))/s^2

# This nested integral makes the computation expensive. Let me instead
# use the swap/Fubini approach as we did for K_c^(II).

# ALTERNATIVE: Use mpmath ODE solver to compute f2, f3, f4 numerically
# at each r, then extract tail amplitudes.

# This is much simpler. Let me use odefun.

print("\n  Strategy: solve f2, f3, f4 ODEs numerically with mpmath")
print("  and extract tail amplitudes via least-squares on cos/sin basis.\n")

R_MAX = mp.mpf(300)
# We'll solve in the variable u_n = r * f_n for each n.
# u_n'' + u_n = Source_n(r)

# Step 1: u1 = sin(r) (exact)
def f1(r):
    if r < mp.mpf('1e-30'):
        return mp.mpf(1)
    return mp.sin(r) / r

def f1p(r):
    if r < mp.mpf('1e-30'):
        return mp.mpf(0)
    return mp.cos(r)/r - mp.sin(r)/r**2

# Step 2: solve u2'' + u2 = -r*(f1')^2
def source2(r):
    fp = f1p(r)
    return -r * fp**2

# Step 3: solve using variation of parameters via accumulating integrals
# u(r) = sin(r)*C(r) - cos(r)*S(r)
# where C(r) = int_0^r cos(s)*source(s) ds
#       S(r) = int_0^r sin(s)*source(s) ds

# For f2:
print("  Computing f2 tail amplitudes (J_c, J_s)...")
Jc = mp.quadosc(lambda s: mp.cos(s) * source2(s), [0, mp.inf], omega=1)
Js = mp.quadosc(lambda s: mp.sin(s) * source2(s), [0, mp.inf], omega=1)
print(f"    J_c (=B2) = {mp.nstr(Jc, 30)}")
print(f"    J_s       = {mp.nstr(Js, 30)}")
print(f"    Expected B2 = {mp.nstr(B2, 30)}")
print(f"    Expected: pi/8 = {mp.nstr(pi/8, 30)}")
print(f"    A2 = -J_s? Let's check: -J_s = {mp.nstr(-Js, 30)}")
# Convention: u2 ~ Jc*sin(r) - Js*cos(r)
# So: sin-amp = Jc = B2, cos-amp = -Js = A2 (if Js = -pi/8... hmm)
# Actually Js should be pi/8, and cos-amp = -Js = -pi/8... but we said A2 = pi/8
# Let me check: I_sin was defined as int sin(s)*(-s*(f1')^2) ds
# source2(s) = -s*(f1')^2, so Js = int sin(s)*source2(s) ds = int sin(s)*(-s*(f1')^2) ds = I_sin
# And I_sin = -pi/8. So Js = I_sin = -pi/8, and A2 = -I_sin = pi/8.
# u2 ~ Jc*sin(r) + (-Js)*(-cos(r))... let me re-derive.

# From var of params: u2(r) = sin(r)*int_0^r cos(s)*src(s)ds - cos(r)*int_0^r sin(s)*src(s)ds
# As r->inf: u2 ~ sin(r)*Jc - cos(r)*Js
# So sin-amplitude = Jc, cos-amplitude = -Js
# u2 = Jc*sin(r) + (-Js)*cos(r) [hmm, let me be precise]
# u2 ~ Jc*sin(r) - Js*cos(r)
# Writing as A*cos + B*sin: B = Jc, A = -Js
# So B2 = Jc, A2 = -Js

print(f"\n    Convention check: B2 = Jc = {mp.nstr(Jc, 20)}, A2 = -Js = {mp.nstr(-Js, 20)}")
A2_val = -Js
B2_val = Jc

# Now for f2 at arbitrary r, we need the full function.
# Instead of computing u2(r) at every point (expensive), let me use
# mpmath odefun to solve the ODEs as a system.

# Actually, let me solve all the ODEs numerically using mpmath's ODE solver.

print("\n  Setting up coupled ODE system...")
print("  Solving u2'' + u2 = src2(r) via mpmath.odefun...")

# mpmath.odefun solves y' = F(r, y) with y(r0) = y0
# For u'' + u = S(r), let w = [u, u']:
# w' = [u', -u + S(r)]

# u2(0) = r*f2(0) = 0 (since f2 is regular at origin)
# u2'(0) = f2(0) + 0*f2'(0) = f2(0)
# From the ODE at r=0: f2(r) ~ f2(0) + f2'(0)*r + f2''(0)*r^2/2 + ...
# f2''(0) = -(f1')^2(0) - 2*f2'(0)/0*... we need L'Hopital for the (2/r)f' term
# At r=0: the ODE f'' + (2/r)f' + f = -(f1')^2
# becomes 3*f''(0) + f(0) = -(f1'(0))^2 = 0 (since f1'(0) = 0)
# So f2''(0) = -f2(0)/3
# For u2 = r*f2: u2(0) = 0, u2'(0) = f2(0)
# What is f2(0)? From the ODE at r=0: 3f2''(0) + f2(0) = 0, so f2''(0) = -f2(0)/3
# The bc is f2(0) = finite, f2'(0) = 0 (regularity)
# Wait, we need to be more careful. The ODE f'' + (2/r)f' + f = -(f1')^2
# with regularity at r=0 determines f2 uniquely (up to adding multiples of f1,
# which we set to zero since the initial condition g0 = 1+delta already
# accounts for the f1 component).
# So f2(0) is determined by the equation. At r=0: f1'(0)=0, so source=0.
# The regular solution has f2'(0)=0, f2(0)=arbitrary (fixed by removing
# homogeneous component). Actually for the perturbation expansion
# g(0) = 1+delta = 1+delta*f1(0)+delta^2*f2(0)+...
# Since g(0) = 1+delta and f1(0) = 1, we need f2(0) = 0.
# Similarly f3(0) = 0, f4(0) = 0.

# So: u2(0) = 0, u2'(0) = f2(0) = 0

# For u2: let me verify by solving with mpmath.
# The source near r=0: source2(r) = -r*(f1')^2
# f1'(r) = cos(r)/r - sin(r)/r^2 ~ -r/3 + ... near r=0
# So source2(r) ~ -r*r^2/9 ~ -r^3/9 near r=0 → source vanishes at r=0.

# Let's solve the system of u2, u3, u4 equations simultaneously.
# But first I need u2 to compute the source for u3, and u3 for u4.
# So let me solve them sequentially.

# Actually, I can solve them all together as a big system:
# State: [u2, u2', u3, u3', u4, u4']

def compute_sources(r, u2, u2p):
    """Compute source terms for the 3rd and 4th order equations."""
    f1_val = f1(r)
    f1p_val = f1p(r)

    # f2 = u2/r, f2' = (u2'*r - u2)/r^2
    if r < mp.mpf('1e-20'):
        f2_val = mp.mpf(0)
        f2p_val = mp.mpf(0)
    else:
        f2_val = u2 / r
        f2p_val = (u2p * r - u2) / r**2

    # Source for f3: -(2*f1'*f2' - f1*(f1')^2)
    src3 = -(2*f1p_val*f2p_val - f1_val*f1p_val**2)
    # For u3: u3'' + u3 = r*src3 (note: fn'' + (2/r)fn' + fn = src3
    # and un = r*fn gives un'' + un = r*src3)
    source3 = r * src3

    return source3, f2_val, f2p_val

# Solve u2 first, then u3, then u4.
# Using mpmath.odefun for each.

# For efficiency, solve all at once as a coupled system.
# But mpmath.odefun may be slow for large r.

# Let me try a direct approach: solve u2 with quadrature (variation of params),
# then u3, then u4, extracting only the tail amplitudes.

# For u2, the tail amplitudes are exactly Jc and Js (already computed).
# For u3, I need: source3(r) = r*[-(2*f1'*f2' - f1*(f1')^2)]
# where f2 = u2(r)/r and f2' = (u2'(r) - u2(r)/r)/r

# This requires u2(r) at every point. I'll compute it via the integral:
# u2(r) = sin(r)*C2(r) - cos(r)*S2(r)
# where C2(r) = int_0^r cos(s)*source2(s) ds
#       S2(r) = int_0^r sin(s)*source2(s) ds

# Then u2'(r) = cos(r)*C2(r) + sin(r)*S2(r)
# (since the differentiation of sin(r)*C2(r) - cos(r)*S2(r) gives
#  cos(r)*C2 + sin(r)*source2*cos - (-sin(r)*S2 + cos(r)*source2*sin)
#  wait, that's not right. Let me redo:
#  d/dr[sin(r)*C2(r)] = cos(r)*C2(r) + sin(r)*cos(r)*source2(r)
#  d/dr[-cos(r)*S2(r)] = sin(r)*S2(r) - cos(r)*sin(r)*source2(r)
#  Sum: u2'(r) = cos(r)*C2(r) + sin(r)*S2(r)
#  ✓ (the source2 terms cancel)

# So: u2(r) = sin(r)*C2(r) - cos(r)*S2(r)
#     u2'(r) = cos(r)*C2(r) + sin(r)*S2(r)

# And f2(r) = u2(r)/r, f2'(r) = (u2'(r)*r - u2(r))/r^2

# For the tail of u3, I need:
# B3 = int_0^inf cos(s)*source3(s) ds  (sin-amplitude = B3)
# A3 = -int_0^inf sin(s)*source3(s) ds  (cos-amplitude = -int sin*src ds)

# Wait, using the same convention:
# u3 ~ sin(r)*int_0^inf cos(s)*src3(s)ds - cos(r)*int_0^inf sin(s)*src3(s)ds
# = sin(r)*B3_check - cos(r)*something

# I need to evaluate:
# int_0^inf cos(s) * source3(s) ds
# int_0^inf sin(s) * source3(s) ds
# where source3(s) = s*[-(2*f1'(s)*f2'(s) - f1(s)*(f1'(s))^2)]

# But this requires f2'(s) at every s, which itself requires C2(s) and S2(s).
# This is the nested/swap integral.

# The swap approach we used for K_c^(II) already handles this!
# K_c^(II) = int_0^inf cos(s)*s*f1'(s)*f2'(s) ds
# which is part of source3.

# source3(s) = s*[-(2*f1'*f2' - f1*(f1')^2)]
# = -2*s*f1'*f2' + s*f1*(f1')^2

# So:
# int cos(s)*source3(s) ds = -2*int cos(s)*s*f1'*f2' ds + int cos(s)*s*f1*(f1')^2 ds
#                          = -2*K_c^(II) + K_c^(I) = P_cos  ✓ (by definition!)

# Similarly:
# int sin(s)*source3(s) ds = -2*int sin(s)*s*f1'*f2' ds + int sin(s)*s*f1*(f1')^2 ds
#                          = -2*K_s^(II) + K_s^(I) = P_sin

# So I need K_s^(I) and K_s^(II):
# K_s^(I) = int_0^inf sin(s)*s*f1(s)*(f1'(s))^2 ds
# K_s^(II) = int_0^inf sin(s)*s*f1'(s)*f2'(s) ds

print("\n  --- Computing P_sin = K_s^(I) - 2*K_s^(II) ---")
print("  Step 1: K_s^(I) = int sin(s)*s*f1*(f1')^2 ds ...")

KsI = mp.quadosc(lambda s: mp.sin(s) * s * f1(s) * f1p(s)**2, [0, mp.inf], omega=1)
print(f"    K_s^(I) = {mp.nstr(KsI, 30)}")

# K_s^(II) requires the swap approach, similar to K_c^(II).
# K_s^(II) = int_0^inf sin(s)*s*f1'(s)*h'(s) ds
# where h = f2, h'(s) = f2'(s) = (u2' - u2/s)/s

# Using the swap: h'(s) = -(cos(s)*Jc(s) + sin(s)*Js(s))/s + (sin(s)*Jc(s) - cos(s)*Js(s))/s^2
# Wait, this requires the partial integrals Jc(s), Js(s).

# Actually, let me use the same Fubini approach as before.
# f2'(r) = u2'(r)/r - u2(r)/r^2
# u2(r) = sin(r)*C2(r) - cos(r)*S2(r)
# u2'(r) = cos(r)*C2(r) + sin(r)*S2(r)
# f2'(r) = [cos(r)*C2(r) + sin(r)*S2(r)]/r - [sin(r)*C2(r) - cos(r)*S2(r)]/r^2
# = [r*cos(r)*C2(r) + r*sin(r)*S2(r) - sin(r)*C2(r) + cos(r)*S2(r)] / r^2

# K_s^(II) = int_0^inf sin(s)*s*f1'(s)*f2'(s) ds
# = int_0^inf sin(s) * f1'(s) * [cos(s)*C2(s)*s + sin(s)*S2(s)*s - sin(s)*C2(s) + cos(s)*S2(s)] / s ds

# This is messy. Let me use the swap:
# C2(s) = int_0^s cos(t)*source2(t) dt
# S2(s) = int_0^s sin(t)*source2(t) dt

# K_s^(II) = int_0^inf sin(s)*s*f1'(s) * [complex expression with C2, S2] / s ds

# Following the same swap trick as for K_c^(II):
# Swap order: K_s^(II) = int_0^inf source2(t) * [int_t^inf kernel(s,t) ds] dt
# where the kernel involves sin(s)*f1'(s)*(cos(t) and sin(t) components)

# This is exactly the same structure as the K_c^(II) computation but with
# sin(s) instead of cos(s) in the outer integral. The tail functions Phi_i
# will be different.

# Let me just use mpmath's ODE solver to compute f2 numerically and then
# do a direct quadosc for K_s^(II).

# Actually, the most straightforward approach: use the ODE solver for the
# full 6-component system [u2, u2', u3, u3', u4, u4'] with mpmath.

print("\n  --- Direct ODE approach with mpmath ---")
print("  Solving u2''+u2=src2, u3''+u3=src3, u4''+u4=src4")
print("  simultaneously from r=eps to r=R_max...")

# Reduce dps for speed
mp.mp.dps = 25
R_MAX = 200
N_steps = 40000
dr = mp.mpf(R_MAX) / N_steps
eps = mp.mpf('1e-8')

# Initialize
r = eps
u2 = mp.mpf(0)
u2p = mp.mpf(0)
u3 = mp.mpf(0)
u3p = mp.mpf(0)
u4 = mp.mpf(0)
u4p = mp.mpf(0)

# RK4 stepper
def rk4_step(r, state, dr):
    u2, u2p, u3, u3p, u4, u4p = state

    def derivs(r, u2, u2p, u3, u3p, u4, u4p):
        f1v = f1(r)
        f1pv = f1p(r)

        # Source for u2
        s2 = -r * f1pv**2

        # f2, f2' from u2
        if r < mp.mpf('1e-15'):
            f2v = mp.mpf(0)
            f2pv = mp.mpf(0)
        else:
            f2v = u2 / r
            f2pv = (u2p * r - u2) / r**2

        # Source for u3: r * [-(2*f1'*f2' - f1*(f1')^2)]
        s3 = r * (-(2*f1pv*f2pv - f1v*f1pv**2))

        # f3, f3' from u3
        if r < mp.mpf('1e-15'):
            f3v = mp.mpf(0)
            f3pv = mp.mpf(0)
        else:
            f3v = u3 / r
            f3pv = (u3p * r - u3) / r**2

        # Source for u4: r * [-((f2')^2 + 2*f1'*f3' - 2*f1*f1'*f2' - f2*(f1')^2 + f1^2*(f1')^2)]
        s4 = r * (-((f2pv)**2 + 2*f1pv*f3pv - 2*f1v*f1pv*f2pv - f2v*f1pv**2 + f1v**2*f1pv**2))

        return [u2p, -u2 + s2, u3p, -u3 + s3, u4p, -u4 + s4]

    k1 = derivs(r, u2, u2p, u3, u3p, u4, u4p)
    s1 = [u2+dr/2*k1[0], u2p+dr/2*k1[1], u3+dr/2*k1[2], u3p+dr/2*k1[3], u4+dr/2*k1[4], u4p+dr/2*k1[5]]
    k2 = derivs(r+dr/2, *s1)
    s2 = [u2+dr/2*k2[0], u2p+dr/2*k2[1], u3+dr/2*k2[2], u3p+dr/2*k2[3], u4+dr/2*k2[4], u4p+dr/2*k2[5]]
    k3 = derivs(r+dr/2, *s2)
    s3 = [u2+dr*k3[0], u2p+dr*k3[1], u3+dr*k3[2], u3p+dr*k3[3], u4+dr*k3[4], u4p+dr*k3[5]]
    k4 = derivs(r+dr, *s3)

    new_state = []
    for i in range(6):
        new_state.append(state[i] + dr/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]))
    return new_state

# Integrate
print(f"  Integrating with RK4, R_max={R_MAX}, steps={N_steps}, dr={mp.nstr(dr,6)}")
state = [u2, u2p, u3, u3p, u4, u4p]

# Collect data for tail extraction
r_data = []
u2_data = []
u3_data = []
u4_data = []

save_interval = max(1, N_steps // 2000)  # save ~2000 points

for step in range(N_steps):
    state = rk4_step(r, state, dr)
    r += dr

    if step % save_interval == 0 and float(r) > 50:
        r_data.append(float(r))
        u2_data.append(float(state[0]))
        u3_data.append(float(state[2]))
        u4_data.append(float(state[4]))

    if step % (N_steps // 5) == 0:
        print(f"    r = {mp.nstr(r, 8)}, u2 = {mp.nstr(state[0], 15)}, u3 = {mp.nstr(state[2], 15)}, u4 = {mp.nstr(state[4], 15)}")

print(f"  Integration complete. R_final = {mp.nstr(r, 10)}")

# Extract tail amplitudes via least-squares
import numpy as np

r_arr = np.array(r_data)
u2_arr = np.array(u2_data)
u3_arr = np.array(u3_data)
u4_arr = np.array(u4_data)

# Use tail region [80, R_MAX-20]
mask = (r_arr > 80) & (r_arr < R_MAX - 20)
r_fit = r_arr[mask]

X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])

def fit_tail(u_arr, name):
    u_fit = u_arr[mask]
    coef, *_ = np.linalg.lstsq(X, u_fit, rcond=None)
    A = coef[0]  # cos amplitude
    B = coef[1]  # sin amplitude
    print(f"    {name}: A (cos) = {A:.15f}, B (sin) = {B:.15f}")
    return A, B

print(f"\n  Tail extraction (fitting window [{r_fit[0]:.0f}, {r_fit[-1]:.0f}]):")
A2_fit, B2_fit = fit_tail(u2_arr, "u2")
A3_fit, B3_fit = fit_tail(u3_arr, "u3")
A4_fit, B4_fit = fit_tail(u4_arr, "u4")

# Check against known values
print(f"\n  Verification against known values:")
print(f"    B2 = {B2_fit:.12f}  (expected I_cos = {float(B2):.12f}, diff = {B2_fit-float(B2):.3e})")
print(f"    A2 = {A2_fit:.12f}  (expected pi/8 = {float(A2):.12f}, diff = {A2_fit-float(A2):.3e})")
print(f"    B3 = {B3_fit:.12f}  (expected P_cos = {float(B3_known):.12f}, diff = {B3_fit-float(B3_known):.3e})")
# A3 = P_sin (new!)
print(f"    A3 = {A3_fit:.12f}  (= -P_sin, NEW)")
print(f"    B4 = {B4_fit:.12f}  (= Q_cos, NEW)")
print(f"    A4 = {A4_fit:.12f}  (= -Q_sin, NEW)")

# Compute alpha_5 from the formula:
# alpha_5 = A3^2/2 + A2*A4 - B2*A2*A3 + A2^2*(B2^2-B3)/2 - A2^4/8
a = B2_fit  # B2
b = B3_fit  # B3
p = A2_fit  # A2 (cos amplitude)
q = A3_fit  # A3
r_val = A4_fit  # A4

# Wait, let me re-check: in my derivation,
# a = B2, b = B3, c = B4, p = A2, q = A3, r = A4
# alpha_5 = q^2/2 + p*r - a*p*q + p^2*(a^2-b)/2 - p^4/8

c = B4_fit
alpha5_formula = q**2/2 + p*r_val - a*p*q + p**2*(a**2-b)/2 - p**4/8
print(f"\n  alpha_5 from formula:")
print(f"    alpha_5 = q^2/2 + p*r - a*p*q + p^2*(a^2-b)/2 - p^4/8")
print(f"    = {alpha5_formula:.10f}")

# Also compute directly: alpha_3 from same amplitudes for verification
alpha3_formula = b + p**2/2
print(f"\n  Verification: alpha_3 = B3 + A2^2/2 = {alpha3_formula:.12f}")
print(f"  Expected alpha_3 = {float(mp.mpf('0.08972222367362532604749')):.12f}")

# alpha_2 from formula:
alpha2_formula = a  # B2 = I_cos
print(f"  alpha_2 = B2 = {alpha2_formula:.12f}")
print(f"  Expected: c1/2 = {float(B2):.12f}")

# alpha_4 from formula:
# d^3 coefficient: c + pq - ap^2/2
alpha4_formula = c + p*q - a*p**2/2
print(f"\n  Bonus: alpha_4 = {alpha4_formula:.10f}")
print(f"  (odd coefficient, from eta_exc - eta_def)")

print(f"\n{'='*78}")
print(f"  FINAL RESULT:")
print(f"    alpha_5 = {alpha5_formula:.10f}")
print(f"    (compare ODE Richardson: ~0.0275)")
print(f"{'='*78}")
