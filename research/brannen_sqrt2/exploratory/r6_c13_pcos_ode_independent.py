#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_c13_pcos_ode_independent.py:
INDEPENDENT verification of alpha_3 = pi^2/128 + P_cos.

Strategy: Bypass the swap entirely. Solve the h-equation as an IVP:
    h'' + (2/r) h' + h = -(f'(r))^2,   h(0)=h_0, h'(0)=0, h bounded at infinity.

Use the closed-form Green's function for the bounded solution:
    h(r) = -u(r)/r,   u(r) = sin(r) J_c(r) - cos(r) J_s(r)
where
    J_c(r) = int_0^r cos(t) t (f'(t))^2 dt
    J_s(r) = int_0^r sin(t) t (f'(t))^2 dt.

Then directly evaluate:
    P_cos = int_0^inf cos(r) r [ f(r)(f'(r))^2 - 2 f'(r) h'(r) ] dr
with h'(r) computed from the same Green's formulas:
    u'(r) = cos(r) J_c(r) + sin(r) J_s(r)
    h'(r) = -u'(r)/r + u(r)/r^2

This avoids the double integral swap entirely — if this matches the swap result to
many digits, the swap path is correct. If it matches pi^2/110 we have closure; if
not, alpha_3 != pi^2/110 exactly.

Method: precompute J_c, J_s on a fine grid via mpmath high-precision cumulative
quadrature, spline-interpolate, then integrate the full P_cos integrand via quadosc.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import mpmath as mp

mp.mp.dps = 35
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

print("="*72)
print(f"  INDEPENDENT alpha_3 CHECK via h-Green formula (dps={mp.mp.dps})")
print("="*72)

def f_fun(s):
    if s < mp.mpf('1e-15'):
        return mp.mpf(1) - s**2/6 + s**4/120 - s**6/5040
    return mp.sin(s)/s

def fp(s):
    if s < mp.mpf('1e-3'):
        return -s/3 + s**3/30 - s**5/840 + s**7/45360 - s**9/3991680
    return (s*mp.cos(s) - mp.sin(s))/s**2

def Aker(t):
    return t * fp(t)**2

# === Step 1: Build J_c(r), J_s(r) via ODE (same precision throughout). ===
# dJ_c/dr = cos(r) r (f')^2,  J_c(0)=0
# dJ_s/dr = sin(r) r (f')^2,  J_s(0)=0
# Use mp.odefun — it gives us an interpolating callable.

print("\n  Building J_c, J_s via mp.odefun ...")

def jc_rhs(r, y): return [mp.cos(r) * Aker(r), mp.sin(r) * Aker(r)]

# odefun accepts vector-valued: y = [J_c, J_s]
# Build up to some large R.
R_max = mp.mpf(50)
print(f"    integrating J_c, J_s on [0, {R_max}] ...", flush=True)
jvec = mp.odefun(jc_rhs, mp.mpf(0), [mp.mpf(0), mp.mpf(0)])
# Probe a value to confirm lazy evaluation works
_probe = jvec(mp.mpf(1))
print(f"    J_c(1), J_s(1) = {mp.nstr(_probe[0], 20)}, {mp.nstr(_probe[1], 20)}")

def Jc(r):
    return jvec(r)[0]
def Js(r):
    return jvec(r)[1]

# === Step 2: h(r) and h'(r) via Green's formula. ===
# Expansions at r->0: J_c(r) ~ r^5/90 (since Aker ~ r^3/9, cos ~ 1),
# J_s(r) ~ r^6/54 (sin ~ r). sin(r)J_c - cos(r)J_s ~ r^6/90 - r^6/54. Small.
# To avoid catastrophic cancellation in u/r at small r, switch to Taylor.
# At r->0, h(r) = h_0 + O(r^2). h_0 = -1/30 (known from previous work).

def u(r):
    return mp.sin(r) * Jc(r) - mp.cos(r) * Js(r)

def up(r):
    return mp.cos(r) * Jc(r) + mp.sin(r) * Js(r)

def h_of_r(r):
    if r < mp.mpf('1e-3'):
        # h(r) near 0: use series expansion from the ODE.
        # h(0) = ? Determined by boundedness. From g expansion, h(0) should be finite.
        # Let me leave computation: for r>=1e-3 the direct formula -u/r is fine.
        r2 = r*r
        # use direct for small too, just watch precision
        pass
    ur = u(r)
    return -ur/r

def hp_of_r(r):
    ur = u(r)
    upr = up(r)
    return -upr/r + ur/r**2

# Quick probe
for rv in [0.01, 0.1, 1.0, 5.0, 20.0]:
    r_mp = mp.mpf(rv)
    hv = h_of_r(r_mp)
    hpv = hp_of_r(r_mp)
    print(f"    r={rv:6.2f}: h = {mp.nstr(hv, 15):>20}, h' = {mp.nstr(hpv, 15):>20}")

# === Step 3: P_cos integrand. ===
def P_integrand(r):
    if r < mp.mpf('1e-8'):
        # integrand ~ r * cos(r) * [f(f')^2 - 2 f' h'] = r * 1 * [1·(r^2/9) - 2·(-r/3)·h'(0)]
        # Leading: r^3/9 + (2r/3) h'(0). But h'(0) = 0 (bounded). So ~ r^3/9, finite.
        return mp.mpf(0)
    fr = f_fun(r)
    fpr = fp(r)
    hpr = hp_of_r(r)
    return mp.cos(r) * r * (fr * fpr**2 - 2 * fpr * hpr)

# === Step 4: Evaluate P_cos via quadosc on [0, R_max] then manage tail. ===
# The ODE solver built J_c, J_s up to R_max=50. Beyond that we need tail extension.
# For r -> inf: J_c, J_s tend to constants (J_c(inf) = ln(3)/8 - 1/2 derived earlier).
# So u(r) -> sin(r) J_c(inf) - cos(r) J_s(inf) (oscillatory bounded),
# h(r) ~ O(1/r), h'(r) ~ O(1/r). f' ~ cos(r)/r.
# Integrand: cos(r) r f' h' ~ cos(r) r · (cos/r) · (osc/r) = cos^2(r) · osc / r, integrable osc.

# We'll compute on [0, R_max] via quadrature on intervals of length pi/2.
print(f"\n  Computing P_cos on [0, {R_max}] ...")
import time
t0 = time.time()

# Split into many sub-intervals for quadrature accuracy
n_subs = 200
edges = [mp.mpf(k)*R_max/n_subs for k in range(n_subs+1)]
Pcos_partial = mp.quad(P_integrand, edges)
print(f"    P_cos on [0, {R_max}] = {mp.nstr(Pcos_partial, 25)}  (took {time.time()-t0:.1f} s)")

# === Step 5: Tail contribution on [R_max, inf]. ===
# Use the asymptotic form of J_c, J_s: they approach constants, so u(r) = sin(r) Jc_inf - cos(r) Js_inf
# plus a rapidly decaying correction from the tail of the Aker·cos/sin integrand.
# At R_max=50, J_c, J_s are essentially at their infinity values modulo ~O(1/R_max^2).

# We use analytic tail with u_tail(r) = sin(r) Jc_inf - cos(r) Js_inf,  h_tail = -u_tail/r.
Jc_inf = ln3/8 - mp.mpf(1)/2
# Js(inf) analytical: int_0^inf sin(r) r (f')^2 dr.
# From previous work: this integral = (pi - 2 Si(2))/8. Let's just evaluate numerically.
# Actually, let's just read J_c(R_max), J_s(R_max) from the ODE.
Jc_Rmax = Jc(R_max)
Js_Rmax = Js(R_max)
print(f"    J_c({R_max})   = {mp.nstr(Jc_Rmax, 20)}")
print(f"    J_c(inf) theo = {mp.nstr(Jc_inf, 20)}")
print(f"    J_s({R_max})   = {mp.nstr(Js_Rmax, 20)}")

# Tail integrand using the values at R_max (first approximation for the tail)
def P_integrand_tail(r):
    if r > mp.mpf('1e8'): return mp.mpf(0)
    # Assume J_c(r) ~ Jc_Rmax, J_s(r) ~ Js_Rmax (frozen) — fine for oscillatory tail
    ur = mp.sin(r) * Jc_Rmax - mp.cos(r) * Js_Rmax
    upr = mp.cos(r) * Jc_Rmax + mp.sin(r) * Js_Rmax
    hpr = -upr/r + ur/r**2
    fr = f_fun(r)
    fpr = fp(r)
    return mp.cos(r) * r * (fr * fpr**2 - 2 * fpr * hpr)

print(f"\n  Tail [R_max, inf] via quadosc ...")
t0 = time.time()
Pcos_tail = mp.quadosc(P_integrand_tail, [R_max, mp.inf], period=2*pi)
print(f"    P_cos tail = {mp.nstr(Pcos_tail, 25)}  (took {time.time()-t0:.1f} s)")

Pcos_total = Pcos_partial + Pcos_tail
print(f"\n  P_cos (total, independent) = {mp.nstr(Pcos_total, 25)}")

alpha3_indep = pi**2/128 + Pcos_total
print(f"  alpha_3 (independent)      = {mp.nstr(alpha3_indep, 25)}")
print(f"  pi^2/110  target           = {mp.nstr(pi**2/110, 25)}")
diff = alpha3_indep - pi**2/110
print(f"  diff                       = {mp.nstr(diff, 10)}")

# Also compare to the swap-based numeric
# From previous work: alpha_3_swap = 0.089722223673625326...
alpha3_swap_prev = mp.mpf('0.089722223673625326047494678804839735')
print(f"  alpha_3 (prev swap)        = {mp.nstr(alpha3_swap_prev, 25)}")
print(f"  diff (indep - swap)        = {mp.nstr(alpha3_indep - alpha3_swap_prev, 10)}")

print("="*72)
print("  Interpretation:")
print("    If (indep - swap) ~ 0 to many digits AND (indep - pi^2/110) ~ -1.5e-6")
print("    => alpha_3 is NOT exactly pi^2/110; the conjecture needs revision.")
print("    If (indep - pi^2/110) ~ 0 to many digits,")
print("    => swap integration has a systematic tail error; revisit quadosc settings.")
print("="*72)
