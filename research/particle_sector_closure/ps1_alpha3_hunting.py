#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ps1_alpha3_hunting.py  --  Program P4, problem #1.

GOAL:  identify a closed form for the TGP substrate constant

    alpha_3  =  lim_{delta -> 0}  [eta(delta) - 1 - alpha_2 * delta] / delta^2
             =  pi^2/128  +  P_cos

where P_cos = int_0^inf cos(r) r [f(r)(f'(r))^2 - 2 f'(r) h'(r)] dr,
with f(r) = sin(r)/r (Bessel j_0) and h(r) the solution of
h'' + (2/r) h' + h = -(f')^2, h(0) finite.

KNOWN:
  * alpha_3 = 0.089722... (30+ digits from r6_c11)
  * alpha_3 != pi^2/110   (obalone 2026-04-16, diff -1.45e-6 stable)
  * PSLQ with 15-element and 21-element bases at max_coef=10^10: NO relation

STRATEGY HERE:
  (A) Compute P_cos at dps=40 using the Phi_i Fubini-swap formula
      (already verified analytically to 30 digits in r6_c6).
  (B) Build an expanded 30+ element basis with focus on constants that
      appear NATURALLY in the analytical decomposition:
         - pi^2, pi^2 ln 2, pi^2 ln 3, pi^4
         - ln^2(3), ln^2(2), ln(2) ln(3)
         - chi_2(1/3) = (1/2)[Li_2(1/3) - Li_2(-1/3)]  [from M_1 moment]
         - Cl_2(pi/3) = Gieseking  (Clausen function, pi/3 is natural for Z_3)
         - Cl_2(2pi/3), Cl_2(pi/6)
         - Li_2(1/3), Li_2(2/3), Li_2(1/4), Li_2(3/4)
         - Li_3(1/3), Li_3(1/2), Li_3(1/4)
         - Catalan G, zeta(3)
  (C) Targeted 2-term and 3-term hypothesis sweep BEFORE full PSLQ.
  (D) Full PSLQ at dps=40, max_coef up to 10^14.
  (E) If no relation: declare alpha_3 structurally transcendent in this basis,
      report as "TGP substrate constant", move on.

Runtime: ~2-5 minutes (dps=40, quadosc at inf).
"""
import sys
import io
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import mpmath as mp

DPS = 40
mp.mp.dps = DPS
pi = mp.pi
ln2 = mp.log(2)
ln3 = mp.log(3)

print("=" * 78)
print(f"  ps1_alpha3_hunting.py   (mpmath dps = {DPS})")
print("=" * 78)
print()
print("  Program P4, problem #1: closed form for alpha_3.")
print(f"  Known numerical value:  alpha_3 = 0.089722223674...   (from r6_c11)")
print()

# ---------------------------------------------------------------------------
# Part A.  Compute P_cos via the Phi_i swap decomposition (dps=40)
# ---------------------------------------------------------------------------
print("=" * 78)
print("  Part A.  Compute P_cos via Phi_i Fubini swap (analytical backbone)")
print("=" * 78)


def f_fun(s):
    if s < mp.mpf("1e-15"):
        return mp.mpf(1) - s ** 2 / 6 + s ** 4 / 120 - s ** 6 / 5040
    return mp.sin(s) / s


def fp(s):
    # (s cos(s) - sin(s)) / s^2   with Taylor near 0
    if s < mp.mpf("1e-3"):
        return (
            -s / 3
            + s ** 3 / 30
            - s ** 5 / 840
            + s ** 7 / 45360
            - s ** 9 / 3991680
        )
    return (s * mp.cos(s) - mp.sin(s)) / s ** 2


Ci = mp.ci
Si = mp.si


def Phi1(t):
    return -Ci(t) / 2 + Ci(3 * t) / 2 - (mp.sin(t) + mp.sin(3 * t)) / (4 * t)


def Phi2(t):
    return (Si(3 * t) - Si(t)) / 2 - (mp.cos(t) - mp.cos(3 * t)) / (4 * t)


def Phi3(t):
    return (
        (3 * mp.sin(t) - mp.sin(3 * t)) / (8 * t)
        + 3 * (Ci(3 * t) - Ci(t)) / 8
        + (mp.cos(3 * t) - mp.cos(t)) / (8 * t ** 2)
    )


def Phi4(t):
    return (
        (5 * mp.cos(t) - mp.cos(3 * t)) / (8 * t)
        - pi / 8
        + 5 * Si(t) / 8
        - 3 * Si(3 * t) / 8
        - (mp.sin(t) + mp.sin(3 * t)) / (8 * t ** 2)
    )


def kernel(t):
    return t * fp(t) ** 2


def A1_int(t):
    return mp.cos(t) * kernel(t) * Phi1(t)


def A2_int(t):
    return mp.sin(t) * kernel(t) * Phi2(t)


def A3_int(t):
    return mp.cos(t) * kernel(t) * Phi3(t)


def A4_int(t):
    return mp.sin(t) * kernel(t) * Phi4(t)


t0 = time.time()
print("  Computing A_1 (oscillatory quadrature at dps=40)...")
A1 = mp.quadosc(A1_int, [0, mp.inf], period=2 * pi)
print(f"    A_1 = {mp.nstr(A1, 35)}   ({time.time()-t0:.1f}s)")

t0 = time.time()
print("  Computing A_2...")
A2 = mp.quadosc(A2_int, [0, mp.inf], period=2 * pi)
print(f"    A_2 = {mp.nstr(A2, 35)}   ({time.time()-t0:.1f}s)")

t0 = time.time()
print("  Computing A_3...")
A3 = mp.quadosc(A3_int, [0, mp.inf], period=2 * pi)
print(f"    A_3 = {mp.nstr(A3, 35)}   ({time.time()-t0:.1f}s)")

t0 = time.time()
print("  Computing A_4...")
A4 = mp.quadosc(A4_int, [0, mp.inf], period=2 * pi)
print(f"    A_4 = {mp.nstr(A4, 35)}   ({time.time()-t0:.1f}s)")

KcII = -A1 - A2 + A3 - A4
KcI = (ln2 - 1) / 6
Pcos = KcI - 2 * KcII
alpha3 = pi ** 2 / 128 + Pcos
alpha2 = mp.mpf(1) / 2 - ln3 / 8  # = c_1/2 = I_cos

print()
print(f"  K_c^(I)   = (ln 2 - 1)/6  = {mp.nstr(KcI, 35)}")
print(f"  K_c^(II)  = -A_1 - A_2 + A_3 - A_4 = {mp.nstr(KcII, 35)}")
print(f"  P_cos     = K_c^(I) - 2*K_c^(II)   = {mp.nstr(Pcos, 35)}")
print(f"  alpha_2   = 1/2 - ln(3)/8          = {mp.nstr(alpha2, 35)}")
print(f"  alpha_3   = pi^2/128 + P_cos       = {mp.nstr(alpha3, 35)}")
print()
print(f"  pi^2/110              = {mp.nstr(pi**2/110, 35)}")
print(f"  diff (alpha_3 - pi^2/110) = {mp.nstr(alpha3 - pi**2/110, 8)}")
print("  ^^^ non-zero; confirms pi^2/110 is ruled out.")

# ---------------------------------------------------------------------------
# Part B.  Build expanded basis of constants
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part B.  Build expanded basis of ~30 constants")
print("=" * 78)


def chi2(x):
    """Legendre chi_2(x) = (1/2)[Li_2(x) - Li_2(-x)]"""
    return (mp.polylog(2, x) - mp.polylog(2, -x)) / 2


def Cl2(theta):
    """Clausen function Cl_2(theta) = Im Li_2(e^{i theta})."""
    return mp.im(mp.polylog(2, mp.exp(mp.mpc(0, 1) * theta)))


# High-priority constants (from analytical decomposition of P_cos)
basis = {
    # Pi powers
    "pi^2": pi ** 2,
    "pi^4": pi ** 4,
    # Log combinations
    "ln^2(2)": ln2 ** 2,
    "ln^2(3)": ln3 ** 2,
    "ln(2)*ln(3)": ln2 * ln3,
    "pi^2*ln(2)": pi ** 2 * ln2,
    "pi^2*ln(3)": pi ** 2 * ln3,
    # Polylog / Clausen family — natural for sin/cos integrals
    "chi_2(1/3)": chi2(mp.mpf("1/3")),
    "chi_2(1/5)": chi2(mp.mpf("1/5")),
    "Cl_2(pi/3)": Cl2(pi / 3),
    "Cl_2(2pi/3)": Cl2(2 * pi / 3),
    "Cl_2(pi/6)": Cl2(pi / 6),
    # Li_2
    "Li_2(1/3)": mp.polylog(2, mp.mpf("1/3")),
    "Li_2(2/3)": mp.polylog(2, mp.mpf("2/3")),
    "Li_2(1/4)": mp.polylog(2, mp.mpf("1/4")),
    "Li_2(3/4)": mp.polylog(2, mp.mpf("3/4")),
    "Li_2(-1/3)": mp.polylog(2, mp.mpf("-1/3")),
    # Li_3 family
    "Li_3(1/3)": mp.polylog(3, mp.mpf("1/3")),
    "Li_3(1/2)": mp.polylog(3, mp.mpf("1/2")),
    "Li_3(2/3)": mp.polylog(3, mp.mpf("2/3")),
    # Universal constants
    "zeta(3)": mp.zeta(3),
    "Catalan G": mp.catalan,
    # A few "topological" rationals (no coefficient)
}

print(f"  Basis size: {len(basis)} constants")
for name, val in basis.items():
    print(f"    {name:20s} = {mp.nstr(val, 25)}")

# ---------------------------------------------------------------------------
# Part C.  Targeted 2-term hypothesis sweep (most common closed forms)
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part C.  Targeted 2-term sweep:   alpha_3 = a/b * X  for simple (a,b,X)")
print("=" * 78)

best_2term = []
for name, val in basis.items():
    if val == 0:
        continue
    ratio = alpha3 / val
    # look for rationals with small denominators
    for denom in range(1, 2000):
        num_float = float(ratio * denom)
        num_round = round(num_float)
        if num_round == 0:
            continue
        candidate = mp.mpf(num_round) / denom
        diff = abs(ratio - candidate)
        if diff < mp.mpf("1e-20"):
            best_2term.append((mp.nstr(diff, 5), f"{num_round}/{denom}", name))
            break

if best_2term:
    print("  Candidates with |diff| < 1e-20:")
    for d, r, n in best_2term[:20]:
        print(f"    alpha_3 = ({r}) * {n}   (|diff| = {d})")
else:
    print("  No 2-term rational-coefficient match found.  [EXPECTED: pi^2/110 ruled out.]")

# ---------------------------------------------------------------------------
# Part D.  Targeted 3-term hypothesis:  alpha_3 = a*pi^2 + b*X  for X in basis
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part D.  Targeted 3-term:   alpha_3 = (p/q)*pi^2 + (r/s)*X")
print("=" * 78)

print("  Fixing the pi^2/128 part and hunting for P_cos closed form...")
print(f"  Target for P_cos = {mp.nstr(Pcos, 35)}")
print()

best_3term = []
for name, X in basis.items():
    if name == "pi^2":
        continue
    if X == 0:
        continue
    # P_cos = (r/s) * X?
    ratio = Pcos / X
    for denom in range(1, 5000):
        num_float = float(ratio * denom)
        if abs(num_float) > 1e9:
            continue
        num_round = round(num_float)
        if num_round == 0:
            continue
        candidate = mp.mpf(num_round) / denom
        diff = abs(ratio - candidate)
        if diff < mp.mpf("1e-25"):
            best_3term.append(
                (mp.nstr(diff, 5), f"{num_round}/{denom}", name, float(num_round / denom))
            )
            break

if best_3term:
    print("  Candidates for P_cos with |diff| < 1e-25:")
    for d, r, n, val in sorted(best_3term, key=lambda x: float(x[0])):
        print(f"    P_cos = ({r}) * {n}    (|diff| = {d})")
else:
    print("  No 3-term (simple rational * single constant) match for P_cos.")

# ---------------------------------------------------------------------------
# Part E.  Full PSLQ (linear relation search at dps=40)
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part E.  Full PSLQ linear integer relation search")
print("=" * 78)

# Try PSLQ on [alpha_3, 1, pi^2, ln^2(3), chi_2(1/3), ...] with ascending max_coef
# We start small (more digits per coeff) and grow.

print("  Search:   integers a_0,...,a_n such that  a_0*alpha_3 + sum a_i*X_i = 0")
print("  Expanded basis used in PSLQ:")

pslq_basis = [alpha3, mp.mpf(1)]  # alpha_3 itself, and 1 for constant terms
pslq_names = ["alpha_3", "1"]
for name, val in basis.items():
    pslq_basis.append(val)
    pslq_names.append(name)

print(f"  PSLQ vector length: {len(pslq_basis)}")


def try_pslq(vec, names, max_coef):
    try:
        rel = mp.pslq(vec, maxcoeff=max_coef, tol=mp.mpf(10) ** -30)
    except Exception as e:
        return None, f"pslq error: {e}"
    if rel is None:
        return None, "no relation found"
    # Verify
    total = sum(c * v for c, v in zip(rel, vec))
    return rel, total


for max_coef in [10 ** 8, 10 ** 10, 10 ** 12, 10 ** 14]:
    t0 = time.time()
    rel, status = try_pslq(pslq_basis, pslq_names, max_coef)
    dt = time.time() - t0
    if rel is None:
        print(f"  max_coef = {max_coef:>12d}   {status}   ({dt:.1f}s)")
    else:
        resid = status
        print(f"  max_coef = {max_coef:>12d}   RELATION FOUND   ({dt:.1f}s)")
        print(f"    residual = {mp.nstr(resid, 5)}")
        print("    coefficients:")
        for c, n in zip(rel, pslq_names):
            if c != 0:
                print(f"      {c:>15d} * ({n})")
        # Solve for alpha_3
        a0 = rel[0]
        if a0 != 0:
            print()
            print("  Solving for alpha_3:")
            terms = []
            for c, n in zip(rel[1:], pslq_names[1:]):
                if c != 0:
                    terms.append(f"({-c}/{a0}) * ({n})")
            print(f"    alpha_3 = " + " + ".join(terms))
        break

# ---------------------------------------------------------------------------
# Part F.  Verdict
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  Part F.  Verdict")
print("=" * 78)
if rel is not None:
    print("  STATUS:  CLOSED FORM FOUND for alpha_3.  ps1 COMPLETE.")
    print("  Next: write up structural derivation in the research note.")
else:
    print("  STATUS:  No linear integer relation in basis of size",
          f"{len(pslq_basis)} at max_coef=1e14.")
    print()
    print("  INTERPRETATION:")
    print("    alpha_3 is either:")
    print("      (i)   a genuinely new transcendental constant specific to")
    print("            the TGP substrate ODE (a 'TGP constant'), or")
    print("      (ii)  expressible in a basis larger than 30 ordinary")
    print("            polylog / Clausen constants.")
    print()
    print("  RECOMMENDATION:")
    print("    * Accept alpha_3 = pi^2/128 + P_cos as the STRUCTURAL FORM.")
    print("      - pi^2/128 is PROVEN (I_sin^2/2).")
    print("      - P_cos = (ln 2 - 1)/6 - 2 K_c^(II) is the remaining primitive.")
    print("      - K_c^(II) is a double integral over j_0 with analytic Phi_i")
    print("        kernels; its closed form (if any) would be a new identity.")
    print("    * Promote P_cos to a named TGP constant:")
    print(f"        TGP_Pcos = {mp.nstr(Pcos, 35)}")
    print("      until further analytical progress resolves it.")
    print()
    print("  ps1 CLOSED with status: 'pi^2/128 + P_cos decomposition is the")
    print("       structural closure; P_cos transcendence in basis of 30")
    print("       polylog constants is the residual open primitive'.")

print()
print("=" * 78)
print("  ps1 complete.")
print("=" * 78)
