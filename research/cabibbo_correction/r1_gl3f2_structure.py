#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
r1_gl3f2_structure.py -- Full algebraic analysis of GL(3, F_2)
==============================================================
Goal: Understand the group structure relevant to Cabibbo correction.

GL(3, F_2) = PSL(2,7) = simple group of order 168.
It acts on F_2^3 (3-dimensional vector space over GF(2)).

Key questions:
  1. What are the conjugacy classes?
  2. How many Z_3 subgroups exist? Are they conjugate?
  3. What is the character table?
  4. How does Z_3 act on the natural 3D representation?
  5. What mixing matrix arises from Z_3 eigenstates?
  6. What sub-leading corrections appear from the full group?

Autor: Claudian (R1 attack)
Data: 2026-04-14
"""

import numpy as np
import galois
from itertools import product
from collections import Counter, defaultdict

# ================================================================
# 1. ENUMERATE GL(3, GF(2))
# ================================================================
print("=" * 80)
print("  GL(3, GF(2)) -- FULL ALGEBRAIC ANALYSIS")
print("=" * 80)

GF2 = galois.GF(2)

def mat_order(M):
    """Order of matrix M in GL(3, GF(2))."""
    I = GF2(np.eye(3, dtype=int))
    power = M.copy()
    for k in range(1, 200):
        if np.array_equal(power, I):
            return k
        power = power @ M  # galois handles mod 2 automatically
    return None

# Generate all 8^3 = 512 possible 3x3 matrices over GF(2),
# keep only invertible ones (det != 0 mod 2)
print("\n[1] Enumerating GL(3, GF(2))...")

def det_gf2(M):
    """Determinant of 3x3 matrix over GF(2) (returns 0 or 1)."""
    m = M.view(np.ndarray).astype(int)
    d = (m[0,0]*(m[1,1]*m[2,2] - m[1,2]*m[2,1])
       - m[0,1]*(m[1,0]*m[2,2] - m[1,2]*m[2,0])
       + m[0,2]*(m[1,0]*m[2,1] - m[1,1]*m[2,0]))
    return d % 2

gl3_elements = []
for bits in range(512):  # 2^9 = 512
    entries = [(bits >> i) & 1 for i in range(9)]
    M = GF2(np.array(entries, dtype=int).reshape(3, 3))
    if det_gf2(M) == 1:
        gl3_elements.append(M)

print(f"  |GL(3, GF(2))| = {len(gl3_elements)}")
assert len(gl3_elements) == 168, f"Expected 168, got {len(gl3_elements)}"

# ================================================================
# 2. ELEMENT ORDERS AND CONJUGACY CLASSES
# ================================================================
print("\n[2] Element orders...")

orders = []
for M in gl3_elements:
    orders.append(mat_order(M))

order_counts = Counter(orders)
print(f"  Order distribution: {dict(sorted(order_counts.items()))}")
# GL(3,F_2) has elements of orders 1, 2, 3, 4, 7

# ================================================================
# 3. CONJUGACY CLASSES
# ================================================================
print("\n[3] Conjugacy classes...")

def mat_to_tuple(M):
    return tuple(M.view(np.ndarray).flatten())

def gf2_inverse(M):
    """Compute inverse of M in GL(3, GF(2)) via adjugate method."""
    m = M.view(np.ndarray).astype(int)
    # Compute adjugate matrix mod 2
    adj = np.zeros((3, 3), dtype=int)
    for i in range(3):
        for j in range(3):
            # Minor (i,j): delete row i, col j
            rows = [r for r in range(3) if r != i]
            cols = [c for c in range(3) if c != j]
            minor = m[np.ix_(rows, cols)]
            cofactor = (minor[0, 0] * minor[1, 1] - minor[0, 1] * minor[1, 0]) % 2
            adj[j, i] = ((-1) ** (i + j) * cofactor) % 2
    # det = 1 (since M is in GL, det=1 mod 2)
    return GF2(adj % 2)

def conjugate(g, h):
    """Compute h^{-1} g h in GF(2)."""
    h_inv = gf2_inverse(h)
    return h_inv @ g @ h

# Build conjugacy classes
assigned = set()
conj_classes = []

for g in gl3_elements:
    g_key = mat_to_tuple(g)
    if g_key in assigned:
        continue
    # Generate conjugacy class of g
    cls = set()
    for h in gl3_elements:
        c = conjugate(g, h)
        cls.add(mat_to_tuple(c))
    for key in cls:
        assigned.add(key)
    # Representative
    rep_order = mat_order(g)
    conj_classes.append({
        'representative': g,
        'order': rep_order,
        'size': len(cls),
        'elements': cls
    })

conj_classes.sort(key=lambda c: (c['order'], c['size']))

print(f"  Number of conjugacy classes: {len(conj_classes)}")
print(f"\n  {'Class':>5s}  {'Order':>5s}  {'Size':>5s}  {'|Centralizer|':>13s}")
print(f"  {'-----':>5s}  {'-----':>5s}  {'-----':>5s}  {'-------------':>13s}")
for i, cls in enumerate(conj_classes):
    centralizer_size = 168 // cls['size']
    print(f"  C{i+1:>3d}    {cls['order']:>3d}    {cls['size']:>3d}    {centralizer_size:>10d}")

# Verify: sum of sizes = 168
total = sum(c['size'] for c in conj_classes)
print(f"\n  Sum of class sizes = {total} (should be 168)")
assert total == 168

# ================================================================
# 4. Z_3 SUBGROUPS
# ================================================================
print("\n[4] Z_3 subgroups (elements of order 3)...")

order3_elements = [M for M, o in zip(gl3_elements, orders) if o == 3]
print(f"  Number of elements of order 3: {len(order3_elements)}")

# Each Z_3 = {I, g, g^2} where g has order 3
# g and g^2 are in the same Z_3, so number of Z_3 subgroups = count/2
z3_subgroups = []
used = set()
I3 = GF2(np.eye(3, dtype=int))

for g in order3_elements:
    g_key = mat_to_tuple(g)
    if g_key in used:
        continue
    g2 = g @ g
    g2_key = mat_to_tuple(g2)
    used.add(g_key)
    used.add(g2_key)
    z3_subgroups.append((g, g2))

print(f"  Number of Z_3 subgroups: {len(z3_subgroups)}")

# Check if all Z_3 subgroups are conjugate
# Two Z_3 = <g> and Z_3' = <g'> are conjugate if g' = h^{-1}gh for some h
z3_conj_classes = []
z3_assigned = set()

for i, (g, g2) in enumerate(z3_subgroups):
    if i in z3_assigned:
        continue
    cls = {i}
    g_key = mat_to_tuple(g)
    g2_key = mat_to_tuple(g2)
    for j, (gp, gp2) in enumerate(z3_subgroups):
        if j in z3_assigned or j == i:
            continue
        gp_key = mat_to_tuple(gp)
        gp2_key = mat_to_tuple(gp2)
        # Check if gp or gp2 is conjugate to g
        for h in gl3_elements:
            c = conjugate(g, h)
            c_key = mat_to_tuple(c)
            if c_key == gp_key or c_key == gp2_key:
                cls.add(j)
                break
    for idx in cls:
        z3_assigned.add(idx)
    z3_conj_classes.append(cls)

print(f"  Z_3 conjugacy classes: {len(z3_conj_classes)}")
for i, cls in enumerate(z3_conj_classes):
    print(f"    Class {i+1}: {len(cls)} Z_3 subgroups")

# ================================================================
# 5. Z_3 EIGENSTATES IN 3D REPRESENTATION
# ================================================================
print("\n[5] Z_3 action on natural 3D representation...")

# Take a representative Z_3 generator
g0 = z3_subgroups[0][0]
print(f"\n  Z_3 generator (over GF(2)):")
print(f"    {g0.view(np.ndarray)}")

# Lift to complex representation
# Over C, a matrix of order 3 has eigenvalues {1, omega, omega^2}
# where omega = e^{2*pi*i/3}
g0_float = g0.view(np.ndarray).astype(float)
omega = np.exp(2j * np.pi / 3)
omega2 = omega**2

eigenvalues, eigenvectors = np.linalg.eig(g0_float)
print(f"\n  Eigenvalues (over R): {eigenvalues}")
print(f"  Note: Over C, order-3 matrix has eigenvalues 1, omega, omega^2")

# Construct the Z_3 action in the generation basis explicitly
# The canonical Z_3 generator acting on 3 generations is the cyclic permutation:
#   P = [[0,0,1],[1,0,0],[0,1,0]]  (cycle 1->2->3->1)
# Its eigenvalues are {1, omega, omega^2}
# Its eigenvectors define the MASS vs FLAVOR basis rotation

P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
print(f"\n  Canonical Z_3 cyclic permutation P:")
print(f"    {P.real.astype(int)}")

eigvals_P, eigvecs_P = np.linalg.eig(P)
print(f"\n  Eigenvalues of P: {eigvals_P}")
print(f"  Expected: 1, omega, omega^2 = 1, {omega:.4f}, {omega2:.4f}")

# Sort eigenvectors by eigenvalue phase
phases = np.angle(eigvals_P)
sort_idx = np.argsort(phases)
eigvals_P = eigvals_P[sort_idx]
eigvecs_P = eigvecs_P[:, sort_idx]

print(f"\n  Eigenvectors (columns = mass eigenstates in flavor basis):")
for i in range(3):
    v = eigvecs_P[:, i]
    print(f"    |m_{i+1}> = [{v[0]:.4f}, {v[1]:.4f}, {v[2]:.4f}]  (eigenvalue {eigvals_P[i]:.4f})")

# ================================================================
# 6. CKM MIXING FROM Z_3 SYMMETRY BREAKING
# ================================================================
print(f"\n{'=' * 80}")
print("  CKM MIXING FROM Z_3 SYMMETRY BREAKING")
print("=" * 80)

# In TGP, the CKM matrix arises from DIFFERENT Z_3 embeddings
# for up-type and down-type quarks.
#
# Key idea: If up-quarks and down-quarks transform under DIFFERENT
# Z_3 subgroups of GL(3,F_2), the CKM matrix is:
#   V_CKM = U_up^dagger @ U_down
# where U_up, U_down diagonalize the respective Z_3 generators.
#
# Two Z_3 generators g_up and g_down in GL(3,F_2) that are NOT conjugate
# (or related by a non-trivial group element) produce a non-trivial CKM.

print("\n  Strategy: V_CKM = U_up^+ @ U_down")
print("  where U_up, U_down diagonalize different Z_3 generators.")

# Take two non-conjugate Z_3 generators (if they exist)
# or two Z_3 generators related by a specific GL(3,F_2) element

# For each pair of Z_3 subgroups, compute the "CKM" mixing
# First, let's parametrize the Z_3 action over C

def z3_to_complex_eigenbasis(g_gf2):
    """
    Given a Z_3 generator over GF(2), lift to a complex unitary matrix
    and return its eigenvector matrix U such that U^+ g U = diag(1, w, w^2).
    """
    g = g_gf2.view(np.ndarray).astype(float)

    # The GF(2) matrix of order 3, lifted to C, has eigenvalues 1, w, w^2
    # But we need to be careful: over R the matrix may not diagonalize
    # Use the formula: for any order-3 matrix M, the projectors are:
    #   P_k = (1/3) sum_{j=0}^{2} omega^{-jk} M^j

    I = np.eye(3)
    M = g.copy()
    M2 = M @ M  # = g^2

    # Projectors onto eigenspaces
    P0 = (I + M + M2) / 3.0          # eigenvalue 1
    P1 = (I + omega**2 * M + omega * M2) / 3.0   # eigenvalue omega
    P2 = (I + omega * M + omega**2 * M2) / 3.0   # eigenvalue omega^2

    # Extract eigenvectors (non-zero columns of projectors)
    vecs = []
    for P_proj, ev in [(P0, 1), (P1, omega), (P2, omega2)]:
        # Find non-zero column
        for col in range(3):
            v = P_proj[:, col]
            if np.abs(v).max() > 1e-10:
                v = v / np.linalg.norm(v)
                vecs.append(v)
                break

    U = np.column_stack(vecs)
    return U

# Compute mixing matrices for all pairs of Z_3 subgroups
print(f"\n  Computing CKM mixing for all {len(z3_subgroups)} x {len(z3_subgroups)} Z_3 pairs...")

mixing_angles = []  # (i, j, theta_12, theta_13, theta_23)

for i, (g_up, _) in enumerate(z3_subgroups):
    U_up = z3_to_complex_eigenbasis(g_up)
    for j, (g_down, _) in enumerate(z3_subgroups):
        if i == j:
            continue
        U_down = z3_to_complex_eigenbasis(g_down)

        # CKM = U_up^+ @ U_down
        V = U_up.conj().T @ U_down

        # Extract |V_us| = |V_{12}| (Cabibbo element)
        V_us = abs(V[0, 1])
        V_ub = abs(V[0, 2])
        V_cb = abs(V[1, 2])

        mixing_angles.append({
            'i': i, 'j': j,
            'V_us': V_us, 'V_ub': V_ub, 'V_cb': V_cb,
            'V': V
        })

# What values of |V_us| appear?
v_us_values = [m['V_us'] for m in mixing_angles]
unique_vus = sorted(set(round(v, 6) for v in v_us_values))
print(f"\n  Distinct |V_us| values from Z_3 pairs: {len(unique_vus)}")
for val in unique_vus:
    count = sum(1 for v in v_us_values if abs(v - val) < 1e-5)
    print(f"    |V_us| = {val:.6f}  ({count} pairs)")

# ================================================================
# 7. THE KEY INSIGHT: RELATIVE ANGLE BETWEEN Z_3 EMBEDDINGS
# ================================================================
print(f"\n{'=' * 80}")
print("  RELATIVE ANGLE BETWEEN Z_3 EMBEDDINGS")
print("=" * 80)

# The relative orientation of two Z_3 subgroups in GL(3,F_2)
# determines the mixing. Over GF(2), the "angle" is discrete.
# But TGP maps this to a continuous mixing via Omega_Lambda.

# Key question: how many DISTINCT relative orientations exist?
# Two pairs (g1, g2) and (g1', g2') give the same CKM if
# there exists h in GL(3,F_2) such that h g1 h^{-1} = g1' and h g2 h^{-1} = g2'

# For now, let's look at the structure constant: for each Z_3 pair,
# what is the "distance" in the group?

print("\n  For each Z_3 pair, computing the group element g_up^{-1} @ g_down...")

# Classify the product g_up^{-1} @ g_down
product_orders = []
for m in mixing_angles:
    i, j = m['i'], m['j']
    g_up = z3_subgroups[i][0]
    g_down = z3_subgroups[j][0]
    g_up_inv = z3_subgroups[i][1]  # g^2 = g^{-1} for order 3
    prod = g_up_inv @ g_down
    prod_ord = mat_order(prod)
    product_orders.append(prod_ord)
    m['product_order'] = prod_ord

prod_order_counts = Counter(product_orders)
print(f"  Order of g_up^{{-1}} @ g_down: {dict(sorted(prod_order_counts.items()))}")

# ================================================================
# 8. PERTURBATIVE CKM FROM SYMMETRY BREAKING
# ================================================================
print(f"\n{'=' * 80}")
print("  PERTURBATIVE CKM FROM Omega_Lambda SYMMETRY BREAKING")
print("=" * 80)

# In TGP, the physical picture is:
# 1. At high energy, GL(3,F_2) is exact -> no mixing (V_CKM = I)
# 2. Omega_Lambda breaks the Z_3 symmetry perturbatively
# 3. The mixing amplitude is proportional to Omega_Lambda/N (tree-level)
# 4. Higher orders give corrections in powers of (Omega_Lambda/N)
#
# The perturbative expansion:
#   V_us = epsilon * <u|Z_3|s> + epsilon^2 * sum_k <u|Z_3|k><k|Z_3|s> / Delta_k + ...
# where epsilon = Omega_Lambda/N
#
# At tree-level: V_us^(1) = Omega_Lambda/N = 0.22823
# At 1-loop:     V_us^(2) = correction from virtual Z_3 transitions
#
# The 1-loop correction involves a sum over intermediate states:
#   delta_V_us = (Omega_Lambda/N)^2 * C_2
# where C_2 is a group-theoretic coefficient.

OMEGA_LAMBDA = 0.6847
N_gen = 3
epsilon = OMEGA_LAMBDA / N_gen

print(f"\n  epsilon = Omega_Lambda/N = {epsilon:.5f}")
print(f"  epsilon^2 = {epsilon**2:.5f}")
print(f"  epsilon^3 = {epsilon**3:.6f}")

# CKM unitarity constraint at O(epsilon^2):
# Sum_k |V_ik|^2 = 1
# |V_ud|^2 + |V_us|^2 + |V_ub|^2 = 1
# At tree-level: V_us = epsilon, V_ub ~ 0 (suppressed), V_ud ~ 1
# -> |V_ud|^2 ~ 1 - epsilon^2
# -> V_ud ~ sqrt(1 - epsilon^2) ~ 1 - epsilon^2/2
#
# But this is the standard Wolfenstein parameterization!
# lambda = sin(theta_C) where:
#   V_ud = cos(theta_C) ~ 1 - lambda^2/2
#   V_us = sin(theta_C) ~ lambda
#
# If TGP gives theta_C = Omega_Lambda/N (the ANGLE, not lambda),
# then lambda = sin(theta_C) = sin(Omega_Lambda/N)
# But sin(0.22823) = 0.22625 -> only 1.9 sigma correction
#
# HOWEVER, there's a deeper point:
# The perturbative expansion of the MIXING AMPLITUDE gives:
#   lambda_C = epsilon - epsilon^3/6 + ... = sin(epsilon)  (trivial)
#
# The non-trivial correction comes from the GROUP STRUCTURE:
#   lambda_C = epsilon * (1 - c * epsilon^2)
# where c is NOT 1/6 (that's just sin) but depends on GL(3,F_2).

# Compute the group-theoretic coefficient c
# From the Z_3 representation theory:
# The number of Z_3 subgroups = n_Z3
# The normalizer N(Z_3) has order |N(Z_3)| = |GL|/n_Z3 * ...

# Let's compute normalizers
print("\n  Computing normalizers of Z_3 subgroups...")

for idx in range(min(3, len(z3_subgroups))):
    g, g2 = z3_subgroups[idx]
    z3_set = {mat_to_tuple(I3), mat_to_tuple(g), mat_to_tuple(g2)}

    normalizer = []
    for h in gl3_elements:
        # h Z_3 h^{-1} = Z_3 iff h g h^{-1} in {g, g^2}
        c = conjugate(g, h)
        c_key = mat_to_tuple(c)
        if c_key in z3_set:
            normalizer.append(h)

    print(f"  Z_3 subgroup #{idx+1}: |Normalizer| = {len(normalizer)}")

    # Normalizer / Z_3 = Aut(Z_3) contribution
    # |Aut(Z_3)| = 2 (phi(3) = 2)
    print(f"    |N(Z_3)/Z_3| = {len(normalizer)//3}")

# ================================================================
# 9. SECOND-ORDER CORRECTION FROM CKM UNITARITY
# ================================================================
print(f"\n{'=' * 80}")
print("  SECOND-ORDER CORRECTION FROM CKM UNITARITY")
print("=" * 80)

# Standard Wolfenstein parameterization to O(lambda^4):
#   V_ud = 1 - lambda^2/2 - lambda^4/8
#   V_us = lambda + O(lambda^7)
#   V_ub = A*lambda^3*(rho - i*eta)
#   V_cd = -lambda + O(lambda^5)
#   V_cs = 1 - lambda^2/2 - lambda^4/8*(1 + 4A^2)
#   V_cb = A*lambda^2
#   V_td = A*lambda^3*(1 - rho - i*eta)
#   V_ts = -A*lambda^2 + O(lambda^4)
#   V_tb = 1 - A^2*lambda^4/2

# In TGP, lambda = Omega_Lambda/N is the leading order.
# The unitarity constraint requires:
#   |V_us|^2 + |V_ud|^2 + |V_ub|^2 = 1
# This is automatically satisfied by the Wolfenstein parameterization.
# But the PHYSICAL question is: does the GROUP STRUCTURE modify lambda?

# Key insight: The Z_3 mixing matrix element is:
#   <u|V|s> = (1/3) * sum_{g in Z_3} omega^{q(g)} * <u|g|s>
# where q(g) is the charge under Z_3 and <u|g|s> is the matrix element.
#
# At tree-level in TGP: <u|g|s> = Omega_Lambda (symmetry breaking parameter)
# At next order: virtual corrections from OTHER group elements
#
# The full sum over GL(3,F_2):
#   V_us = (Omega_Lambda/N) * [1 + sum_{g not in Z_3} c_g * (Omega_Lambda/N)^{n_g}]
#
# The leading correction comes from elements of ORDER 2 (involutions).

# Count elements by order and their contribution
print("\n  Group structure analysis:")
print(f"  Elements of order 1: {order_counts.get(1, 0)} (identity)")
print(f"  Elements of order 2: {order_counts.get(2, 0)} (involutions)")
print(f"  Elements of order 3: {order_counts.get(3, 0)} (Z_3 generators)")
print(f"  Elements of order 4: {order_counts.get(4, 0)} (Z_4 elements)")
print(f"  Elements of order 7: {order_counts.get(7, 0)} (Sylow-7)")

# The correction from group theory:
# At 2nd order, the amplitude gets a correction from the "propagator"
# through intermediate states. In the group-theoretic language:
#
# delta_lambda = - epsilon^2 * (1/|G|) * sum_{g not in Z_3} Tr(g) / (1 - eigenvalue(g))
#
# But this is getting complicated. Let's try a different approach.

# ================================================================
# 10. DIRECT APPROACH: EFFECTIVE N FROM GL(3,F_2) REPRESENTATION THEORY
# ================================================================
print(f"\n{'=' * 80}")
print("  EFFECTIVE N FROM REPRESENTATION THEORY")
print("=" * 80)

# GL(3, F_2) ~ PSL(2,7) has 6 conjugacy classes -> 6 irreps
# Character table of PSL(2,7):
#
# Class:    1A    2A    3A    4A    7A    7B
# Size:      1    21    56    42    24    24
# Order:     1     2     3     4     7     7
#
# chi_1:     1     1     1     1     1     1    (trivial)
# chi_3:     3    -1     0    1    zeta  zeta*  (natural)
# chi_3':    3    -1     0    1   zeta*  zeta   (conjugate)
# chi_6:     6     2     0     0    -1    -1
# chi_7:     7    -1     1    -1     0     0
# chi_8:     8     0    -1     0     1     1
#
# where zeta = (1+i*sqrt(7))/2, zeta* = (1-i*sqrt(7))/2

print("\n  Character table of GL(3,F_2) ~ PSL(2,7):")
print("  (6 conjugacy classes, 6 irreps)")

# Verify our conjugacy classes match
print(f"\n  Our conjugacy classes:")
for i, cls in enumerate(conj_classes):
    print(f"    C{i+1}: order {cls['order']}, size {cls['size']}")

# The NATURAL 3-dimensional representation chi_3 is key.
# Under Z_3 (class 3A), chi_3 decomposes as:
#   chi_3|_{Z_3} = 1 + omega + omega^2 = rho_reg (regular rep of Z_3)
# This means all three Z_3 eigenspaces are 1-dimensional -> good!

# Now, the key formula for the Cabibbo angle correction:
# In the natural rep, the Z_3 action mixes generations.
# The mixing amplitude at tree-level is epsilon = Omega_Lambda/N.
#
# The 1-loop correction involves a sum over virtual states
# weighted by group characters:
#
# V_us = epsilon * [1 + (epsilon^2) * C_2 + (epsilon^4) * C_4 + ...]
#
# where C_2 = -(1/|G|) * sum_{class C} |C| * chi_3(C) * f(C)
#
# The function f(C) encodes the "propagator" in group space.

# For the Z_3 mediated mixing, the relevant correction is:
# C_2 = Tr[sum over non-Z_3 paths] / Tr[Z_3 path]

# Let's compute it directly using the natural 3D rep
print("\n  Computing correction coefficient C_2 from character theory...")

# Characters of chi_3 on each class
chi3_values = {}
for cls in conj_classes:
    rep = cls['representative']
    # Trace in the natural 3D representation over GF(2)
    # Lifted to integers
    tr = int(np.trace(rep.view(np.ndarray).astype(int))) % 2
    # But we need the trace over C, not GF(2)
    # For a GF(2) matrix, the trace over C equals the number of
    # eigenvalue-1 eigenvalues (counted with multiplicity)
    # Use the actual trace of the integer matrix
    tr_int = np.trace(rep.view(np.ndarray).astype(int))
    chi3_values[cls['order']] = chi3_values.get(cls['order'], [])
    chi3_values[cls['order']].append((cls['size'], int(tr_int)))

print(f"  Traces in natural 3D representation:")
for cls in conj_classes:
    rep = cls['representative']
    tr = int(np.trace(rep.view(np.ndarray).astype(int)))
    print(f"    Order {cls['order']}, size {cls['size']}: Tr = {tr}")

# ================================================================
# 11. N_eff FROM GROUP AVERAGING
# ================================================================
print(f"\n{'=' * 80}")
print("  N_eff FROM GROUP AVERAGING")
print("=" * 80)

# Here's a key insight: The effective number of generations might be
# modified by the group structure. Instead of N = 3 (dimension of rep),
# we might need N_eff = <chi_3> averaged appropriately.
#
# Several candidates:
#
# (a) N_eff = (1/|G|) * sum_g |chi_3(g)|^2 = <|chi_3|^2>
#     By orthogonality: sum |chi_3(g)|^2 = |G| = 168  (irreducible)
#     So <|chi_3|^2> = 168/168 = 1 ... not useful
#
# (b) N_eff = (1/|G|) * sum_g chi_3(g) = multiplicity of trivial rep in chi_3
#     = 0 (chi_3 is irreducible and non-trivial)
#
# (c) N_eff = dim(chi_3) * [1 - correction from Casimir]
#     The quadratic Casimir for chi_3: C_2(3) = ...

# More promising: Consider the Z_3 INDEX in GL(3,F_2)
# [GL(3,F_2) : N_G(Z_3)] = 168 / |N_G(Z_3)|
# This counts the number of Z_3 cosets -> related to mixing multiplicity

# From step 4, we found the normalizer size. Let's compute it properly.
g_z3, g2_z3 = z3_subgroups[0]
z3_set_0 = {mat_to_tuple(I3), mat_to_tuple(g_z3), mat_to_tuple(g2_z3)}

normalizer_0 = []
for h in gl3_elements:
    c = conjugate(g_z3, h)
    c_key = mat_to_tuple(c)
    if c_key in z3_set_0:
        normalizer_0.append(h)

N_norm = len(normalizer_0)
n_cosets = 168 // N_norm

print(f"\n  |N_G(Z_3)| = {N_norm}")
print(f"  [G : N_G(Z_3)] = {n_cosets} (number of distinct Z_3 cosets)")
print(f"  Number of Z_3 subgroups = {len(z3_subgroups)}")

# The effective mixing involves N_gen = 3 generations, but the
# number of "mixing channels" is n_cosets.
# Each coset contributes a mixing path with amplitude ~ epsilon.
# But paths interfere, so the total amplitude is:
#   V_us = epsilon * sum of phase factors from cosets

# ================================================================
# 12. CRITICAL CALCULATION: EXACT CKM FROM GL(3,F_2) STRUCTURE
# ================================================================
print(f"\n{'=' * 80}")
print("  CRITICAL: EXACT V_us FROM GL(3,F_2) AVERAGING")
print("=" * 80)

# Let's try the most rigorous approach:
#
# In TGP, the CKM matrix element V_us comes from the amplitude
# for an up-quark (generation 1, Z_3 eigenvalue 1) to transition
# to a strange-quark (generation 2, Z_3 eigenvalue omega).
#
# The amplitude is:
#   A(u -> s) = (epsilon/|G|) * sum_{g in G} chi_3(g) * omega^{-q_3(g)}
#             + higher order
#
# where q_3(g) is the Z_3 charge of g (i.e., the coset index in G/Z_3,
# but this only makes sense for g in N_G(Z_3)).
#
# Actually, the cleaner formulation is:
# The mixing amplitude between generation i and generation j is
# proportional to the (i,j) matrix element of the GROUP AVERAGE:
#
#   M_{ij} = (1/|G|) * sum_{g in G} rho(g)_{ij}
#
# where rho is the 3D representation.
#
# But for irreducible rho: sum_g rho(g) = 0 (Schur's lemma)!
# So the leading term vanishes -- which means the mixing is
# a SYMMETRY-BREAKING effect, as expected.
#
# The symmetry-breaking mixing at order epsilon:
#   V_us = epsilon * <1| delta_H |2>
# where delta_H is the symmetry-breaking Hamiltonian.
#
# In TGP, delta_H is generated by Z_3 breaking.
# The Z_3 breaking lifts the degeneracy of the 3 generations.

# Let's try a concrete model:
# H = H_0 + epsilon * H_1
# H_0 = identity (degenerate)
# H_1 = Z_3 generator (breaks degeneracy)
# Then V_us comes from the off-diagonal elements of the Z_3 generator
# in the PERTURBED basis.

# Z_3 generator in the generation basis:
# P = [[0,0,1],[1,0,0],[0,1,0]]
# Eigenvalues: 1, omega, omega^2
# Eigenvectors: |v_k> = (1/sqrt(3)) * [1, omega^k, omega^{2k}]

# Mass basis -> flavor basis rotation:
U_mass_to_flavor = (1/np.sqrt(3)) * np.array([
    [1, 1, 1],
    [1, omega, omega2],
    [1, omega2, omega]
], dtype=complex)

print("\n  Mass-to-flavor rotation U (from Z_3 diagonalization):")
for i in range(3):
    row = U_mass_to_flavor[i]
    print(f"    [{row[0]:.4f}, {row[1]:.4f}, {row[2]:.4f}]")

# In TGP, up-type quarks and down-type quarks have DIFFERENT
# Z_3 generators (different embeddings in GL(3,F_2)).
# The relative angle between the two Z_3's determines V_CKM.
#
# If the two Z_3's are related by a group element h:
#   Z_3^{down} = h Z_3^{up} h^{-1}
# Then:
#   V_CKM = U_up^+ @ h_rep @ U_down
# where h_rep is the 3D representation of h.
#
# The key: h is NOT arbitrary -- it's a specific element of GL(3,F_2)
# determined by the dynamics (Omega_Lambda breaking).

# Let's enumerate all possible V_CKM matrices
print("\n  All possible |V_us| values from GL(3,F_2) embeddings:")
print("  (considering all distinct h that conjugate one Z_3 to another)")

# For a fixed Z_3^up, enumerate all g such that g Z_3 g^{-1} != Z_3
# (i.e., g NOT in normalizer)
# For each such g, V_CKM = U^+ @ rho(g) @ U

z3_gen = z3_subgroups[0][0]
U_z3 = z3_to_complex_eigenbasis(z3_gen)

vus_from_h = []

for h in gl3_elements:
    h_key = mat_to_tuple(h)
    # Skip identity
    if np.array_equal(h, I3):
        continue

    # Representation of h over C (just the integer matrix)
    h_rep = h.view(np.ndarray).astype(complex)

    # CKM-like matrix: U^+ h U
    V = U_z3.conj().T @ h_rep @ U_z3

    vus = abs(V[0, 1])
    vus_from_h.append((vus, mat_order(h), h_key))

# Distinct values
vus_rounded = Counter(round(v[0], 6) for v in vus_from_h)
print(f"\n  Distinct |V_us| values (from U^+ h U for all h in GL(3,F_2)):")
for val in sorted(vus_rounded.keys()):
    print(f"    |V_us| = {val:.6f}  (count: {vus_rounded[val]})")

# ================================================================
# 13. THE CABIBBO FORMULA WITH GROUP-THEORETIC COEFFICIENT
# ================================================================
print(f"\n{'=' * 80}")
print("  CABIBBO FORMULA WITH GROUP-THEORETIC COEFFICIENTS")
print("=" * 80)

# The TGP formula lambda_C = Omega_Lambda / N is a tree-level result.
# The full result should be:
#   lambda_C = (Omega_Lambda / N) * F(GL(3,F_2))
# where F is a form factor from the group structure.
#
# Candidate form factors:

F_candidates = {
    "1 (tree-level)": 1.0,
    "(|G|-1)/|G| = 167/168": 167/168,
    "(|G|-|Z_3|)/|G| = 165/168": 165/168,
    "sqrt(|G|-1)/sqrt(|G|)": np.sqrt(167/168),
    "1 - 1/|Z_3 orbits|": 1.0 - 1/n_cosets,
    "N/(N+N_eff_corr)": 3.0 / (3.0 + 0.043),  # from ex274
    "|C_3A|/|G| renormalization": 1.0 - 56/168,  # 56 = size of order-3 class
}

LAMBDA_PDG = 0.22500
SIGMA_PDG = 0.00067

print(f"\n  lambda_C = (Omega_Lambda/N) * F")
print(f"  Omega_Lambda/N = {epsilon:.5f}")
print(f"  PDG: {LAMBDA_PDG:.5f} +/- {SIGMA_PDG:.5f}")
print(f"\n  {'Form factor F':<35s}  {'F value':>8s}  {'lambda_C':>10s}  {'Tension':>8s}")
print(f"  {'-'*35}  {'-'*8}  {'-'*10}  {'-'*8}")

for name, F in F_candidates.items():
    lam = epsilon * F
    t = abs(lam - LAMBDA_PDG) / SIGMA_PDG
    star = " <--" if t < 2.0 else ""
    print(f"  {name:<35s}  {F:8.5f}  {lam:10.5f}  {t:6.1f}s{star}")

# ================================================================
# 14. THE MOST PROMISING APPROACH: CKM FROM DOUBLE COSET
# ================================================================
print(f"\n{'=' * 80}")
print("  DOUBLE COSET DECOMPOSITION")
print("=" * 80)

# The CKM matrix is determined by the double coset Z_3 \ GL(3,F_2) / Z_3.
# Each double coset Z_3 g Z_3 gives a distinct mixing pattern.

print("\n  Computing Z_3 \\ GL(3,F_2) / Z_3 double cosets...")

z3_elements_set = {mat_to_tuple(I3), mat_to_tuple(g_z3), mat_to_tuple(g2_z3)}

# Helper: multiply two GF2 matrices
def gf2_mul(A, B):
    return A @ B

# Generate double cosets
assigned_dc = set()
double_cosets = []

for g in gl3_elements:
    g_key = mat_to_tuple(g)
    if g_key in assigned_dc:
        continue

    # Generate Z_3 g Z_3 = {z1 @ g @ z2 : z1, z2 in Z_3}
    dc = set()
    z3_mats = [I3, g_z3, g2_z3]
    for z1 in z3_mats:
        for z2 in z3_mats:
            elem = gf2_mul(gf2_mul(z1, g), z2)
            dc.add(mat_to_tuple(elem))

    for key in dc:
        assigned_dc.add(key)

    double_cosets.append({
        'representative': g,
        'size': len(dc),
        'order': mat_order(g)
    })

double_cosets.sort(key=lambda d: (d['order'], d['size']))

print(f"\n  Number of double cosets: {len(double_cosets)}")
print(f"\n  {'DC':>4s}  {'|DC|':>5s}  {'Order(rep)':>10s}")
print(f"  {'----':>4s}  {'-----':>5s}  {'----------':>10s}")
total_dc = 0
for i, dc in enumerate(double_cosets):
    print(f"  DC{i+1:>2d}  {dc['size']:>5d}  {dc['order']:>10d}")
    total_dc += dc['size']
print(f"\n  Sum of DC sizes = {total_dc} (should be 168)")

# For each double coset, compute the CKM mixing
print(f"\n  CKM mixing from each double coset:")
print(f"  {'DC':>4s}  {'|V_us|':>8s}  {'|V_ub|':>8s}  {'|V_cb|':>8s}  {'Order':>5s}")
print(f"  {'----':>4s}  {'--------':>8s}  {'--------':>8s}  {'--------':>8s}  {'-----':>5s}")

for i, dc in enumerate(double_cosets):
    g = dc['representative']
    g_rep = g.view(np.ndarray).astype(complex)
    V = U_z3.conj().T @ g_rep @ U_z3
    v_us = abs(V[0, 1])
    v_ub = abs(V[0, 2])
    v_cb = abs(V[1, 2])
    dc['V_us'] = v_us
    dc['V_ub'] = v_ub
    dc['V_cb'] = v_cb
    print(f"  DC{i+1:>2d}  {v_us:8.5f}  {v_ub:8.5f}  {v_cb:8.5f}  {dc['order']:>5d}")

# ================================================================
# 15. PHYSICAL INTERPRETATION OF THE DOUBLE COSET STRUCTURE
# ================================================================
print(f"\n{'=' * 80}")
print("  PHYSICAL INTERPRETATION")
print("=" * 80)

print("""
  KEY FINDING: The CKM mixing is determined by which DOUBLE COSET
  of Z_3 in GL(3,F_2) is selected by the dynamics.

  In TGP, the symmetry-breaking parameter Omega_Lambda selects
  a specific double coset (or mixture of cosets weighted by Omega_Lambda).

  The tree-level formula lambda_C = Omega_Lambda/N = 0.22823 corresponds
  to a SPECIFIC group-theoretic structure. The correction comes from
  the WEIGHT of each double coset in the physical amplitude.

  The physical amplitude is:
    V_us = sum_DC w(DC) * V_us(DC)
  where w(DC) is the weight proportional to Omega_Lambda^{n(DC)}.
""")

# ================================================================
# 16. WEIGHTED AVERAGE APPROACH
# ================================================================
print(f"\n{'=' * 80}")
print("  WEIGHTED AVERAGE: V_us FROM GROUP AVERAGING")
print("=" * 80)

# One natural approach: weight by double coset size
total_size = sum(dc['size'] for dc in double_cosets)
V_us_avg_size = sum(dc['size'] * dc['V_us'] for dc in double_cosets) / total_size
V_us_avg_unweighted = np.mean([dc['V_us'] for dc in double_cosets])

print(f"\n  Size-weighted average |V_us|: {V_us_avg_size:.6f}")
print(f"  Unweighted average |V_us|:    {V_us_avg_unweighted:.6f}")
print(f"  PDG lambda:                    {LAMBDA_PDG:.6f}")

# Another approach: The physical lambda_C is:
# lambda_C = (Omega_Lambda/N) * <V_us>_group
# where <V_us>_group is the group-theoretic form factor

# Using each DC's V_us as the form factor:
print(f"\n  lambda_C = (Omega_Lambda/N) * F_DC for each double coset:")
for i, dc in enumerate(double_cosets):
    if dc['V_us'] > 0:
        F = dc['V_us']
        lam = epsilon * F
        # But this doesn't make sense dimensionally...
        # V_us from DC is already a number, not a form factor
        pass

# The CORRECT interpretation:
# TGP predicts lambda_C = Omega_Lambda/N at tree level.
# The GROUP provides a multiplicative correction:
# lambda_C = (Omega_Lambda/N) * [1 - delta_group]
#
# The most natural delta_group from the group structure:

# (A) Identity fraction: The identity DC contributes V_us = 0 (no mixing).
#     Removing it: delta = 1/|G| = 1/168
id_dc = [dc for dc in double_cosets if dc['order'] == 1]
if id_dc:
    print(f"\n  (A) Identity double coset: size = {id_dc[0]['size']}, V_us = {id_dc[0]['V_us']:.6f}")

# (B) Z_3 elements themselves don't mix (diagonal in mass basis).
z3_dc = [dc for dc in double_cosets if dc['order'] == 3]
print(f"  (B) Z_3 double cosets: {len(z3_dc)}")
for dc in z3_dc:
    print(f"      size = {dc['size']}, V_us = {dc['V_us']:.6f}")

# (C) Non-mixing elements fraction:
non_mixing = sum(dc['size'] for dc in double_cosets if dc['V_us'] < 0.01)
mixing = sum(dc['size'] for dc in double_cosets if dc['V_us'] >= 0.01)
print(f"\n  (C) Non-mixing DC elements: {non_mixing}/{total_size}")
print(f"      Mixing DC elements:     {mixing}/{total_size}")
print(f"      Mixing fraction:         {mixing/total_size:.5f}")

correction_factor = mixing / total_size
lambda_corrected = epsilon * correction_factor
print(f"\n  lambda_C (corrected) = epsilon * (mixing/total)")
print(f"                       = {epsilon:.5f} * {correction_factor:.5f}")
print(f"                       = {lambda_corrected:.5f}")
t_corr = abs(lambda_corrected - LAMBDA_PDG) / SIGMA_PDG
print(f"  Tension: {t_corr:.1f} sigma")

# ================================================================
# 17. SUMMARY OF ALGEBRAIC RESULTS
# ================================================================
print(f"\n{'=' * 80}")
print("  SUMMARY OF ALGEBRAIC RESULTS")
print("=" * 80)

print(f"""
  GL(3,F_2) STRUCTURE:
    |G| = 168
    Conjugacy classes: {len(conj_classes)}
    Elements by order: {dict(sorted(order_counts.items()))}

  Z_3 SUBGROUPS:
    Number: {len(z3_subgroups)}
    Conjugacy classes of Z_3: {len(z3_conj_classes)}
    |Normalizer(Z_3)| = {N_norm}
    Z_3 cosets: {n_cosets}

  DOUBLE COSETS Z_3 \\ G / Z_3:
    Number: {len(double_cosets)}
    Sizes: {[dc['size'] for dc in double_cosets]}

  CKM MIXING VALUES (from double cosets):
""")

for i, dc in enumerate(double_cosets):
    print(f"    DC{i+1}: |V_us| = {dc['V_us']:.5f}, order = {dc['order']}, size = {dc['size']}")

print(f"""
  CORRECTION CANDIDATES:
    Tree-level: lambda = Omega_Lambda/N = {epsilon:.5f} (4.8 sigma)
    167/168:    lambda = {epsilon * 167/168:.5f} ({abs(epsilon * 167/168 - LAMBDA_PDG)/SIGMA_PDG:.1f} sigma)
    165/168:    lambda = {epsilon * 165/168:.5f} ({abs(epsilon * 165/168 - LAMBDA_PDG)/SIGMA_PDG:.1f} sigma)
    Mixing frac: lambda = {lambda_corrected:.5f} ({t_corr:.1f} sigma)

  PDG: lambda = {LAMBDA_PDG} +/- {SIGMA_PDG}
""")

print("=" * 80)
print("  END OF GL(3,F_2) ALGEBRAIC ANALYSIS")
print("=" * 80)
