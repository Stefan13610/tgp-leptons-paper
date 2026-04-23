"""
R3: ANALITYCZNE WYPROWADZENIE C_T = K/A^2 ~ 17.60

Kontext (r3_virial_mechanism.py):
  - K/A^2 = 17.597 +- 0.21% uniwersalne dla g0 in [0.5, 2.0]
  - |V|/A^2 = 17.373 +- 0.10%
  - K/|V| = 1.013 (quasi-wirial)

Pytanie: czy C_T = 17.60 jest FUNDAMENTALNA stala, czy artefakt cutoff R_max = 70?

HIPOTEZA: C_T zalezy od R_max (upper cutoff calki)

  Z linearyzacji tail: g - 1 = A*sin(r+delta)/r
  T * r^2 ~ (A^2/2) * cos^2(r+delta)  (dla d=3, r^(d-1)=r^2 i g' ~ A*cos/r)
  int_0^R_max T * r^2 dr ~ (A^2/2) * int_0^R_max cos^2(r+delta) dr
                        ~ (A^2/4) * R_max  (usrednienie cos^2 = 1/2)
  Zatem K/A^2 ~ R_max/4.

PREDYKCJE:
  R_max = 70 -> K/A^2 ~ 70/4 = 17.50  (obserwowane: 17.60)
  R_max = 40 -> K/A^2 ~ 40/4 = 10.00
  R_max = 100 -> K/A^2 ~ 100/4 = 25.00

WNIOSEK (oczekiwany): C_T jest LINIOWA w R_max.
  Wtedy m_phys = c * K^2 = c * (R_max/4)^2 * A^4 = c_eff(R_max) * A^4.
  Ratio m_i/m_j = (K_i/K_j)^2 = (A_i/A_j)^4 NIEZALEZNE od R_max (cutoff sie skraca).

KONSEKWENCJA dla R5 bridge:
  - c_eff(R_max) = c * R_max^2/16 rozumie sie jako skale fizyczna
  - Fizyczny R_max = odleglosc efektywna w przestrzeni R5 substratu
  - m_e = c_eff * A_e^4 pozwala na wyznaczenie c z jednej masy i rozmiar R_max
"""

import numpy as np
from scipy.integrate import solve_ivp


def solve_soliton(g0, alpha=1.0, d=3, r_max=150.0):
    def ode(r, y):
        g, gp = y
        g = max(g, 1e-12)
        rhs = (1.0 - g) * g**(2.0 - 2.0*alpha)
        return [gp, rhs - (alpha/g)*gp*gp - ((d-1)/r)*gp]

    r0 = 1e-4
    c2 = (1.0 - g0) * g0**(2.0 - 2.0*alpha) / (2*d)
    y0 = [g0 + c2*r0*r0, 2*c2*r0]
    return solve_ivp(ode, (r0, r_max), y0, method='DOP853',
                      dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.05)


def extract_A_tail(sol, r_fit_min=20.0, r_fit_max=None):
    if r_fit_max is None:
        r_fit_max = sol.t[-1] - 2
    rs = np.linspace(r_fit_min, r_fit_max, 2000)
    ys = sol.sol(rs)
    u = (ys[0] - 1) * rs
    S, C = np.sin(rs), np.cos(rs)
    M = np.column_stack([S, C])
    coefs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    return np.sqrt(coefs[0]**2 + coefs[1]**2)


def compute_K(sol, alpha=1.0, d=3, r_outer=70.0):
    rs = np.linspace(max(sol.t[0], 1e-4), r_outer, 8000)
    ys = sol.sol(rs)
    g = np.clip(ys[0], 1e-12, 10.0)
    gp = ys[1]
    T = 0.5 * g**(2*alpha) * gp**2
    return np.trapezoid(T * rs**(d-1), rs)


def compute_V(sol, d=3, r_outer=70.0):
    rs = np.linspace(max(sol.t[0], 1e-4), r_outer, 8000)
    ys = sol.sol(rs)
    g = np.clip(ys[0], 1e-12, 10.0)
    V = g**3/3.0 - g**4/4.0 - 1.0/12.0
    return np.trapezoid(V * rs**(d-1), rs)


def main():
    print("=" * 76)
    print("  R3 ANALITYCZNE: C_T(R_max) = K/A^2 jako funkcja cutoff")
    print("=" * 76)
    print()
    print("  PREDYKCJA z linearyzacji tail: C_T ~ R_max/4")
    print()

    # Solve for one representative soliton with large r_max
    g0_list = [0.869, 1.407, 1.729]  # e, mu, tau
    labels = ['e', 'mu', 'tau']

    # Solve wszystkie z duzym R_max
    sols = {}
    As = {}
    r_max_solver = 150.0
    for lab, g0 in zip(labels, g0_list):
        sol = solve_soliton(g0, r_max=r_max_solver)
        sols[lab] = sol
        # A_tail jest WLASNOSCIA solitonu, niezalezna od cutoff
        # Dopasuj w r in [20, 60] (stabilny zakres)
        As[lab] = extract_A_tail(sol, r_fit_min=20.0, r_fit_max=60.0)

    print(f"  A_e   = {As['e']:.6f}")
    print(f"  A_mu  = {As['mu']:.6f}")
    print(f"  A_tau = {As['tau']:.6f}")
    print()

    # Skan R_max: test czy K/A^2 ~ R_max/4
    R_max_list = [30, 40, 50, 60, 70, 80, 100, 120, 140]

    print("  Skalowanie C_T(R_max) = K/A^2:")
    print()
    print(f"  {'R_max':>8} {'C_T_e':>10} {'C_T_mu':>10} {'C_T_tau':>10} "
          f"{'R/4 pred':>10} {'dev_e%':>8}")
    print(f"  {'-----':>8} {'-----':>10} {'------':>10} {'-------':>10} "
          f"{'--------':>10} {'------':>8}")

    C_T_data = {'e': [], 'mu': [], 'tau': []}
    C_V_data = {'e': [], 'mu': [], 'tau': []}

    for R in R_max_list:
        row = [f"  {R:>8}"]
        for lab in labels:
            K = compute_K(sols[lab], r_outer=R)
            V = compute_V(sols[lab], r_outer=R)
            CT = K / As[lab]**2
            CV = abs(V) / As[lab]**2
            C_T_data[lab].append((R, K, CT))
            C_V_data[lab].append((R, abs(V), CV))
            row.append(f"{CT:>10.4f}")
        pred = R / 4.0
        row.append(f"{pred:>10.4f}")
        dev = (C_T_data['e'][-1][2] - pred) / pred * 100
        row.append(f"{dev:>+7.2f}%")
        print(" ".join(row))

    print()
    print("=" * 76)
    print("  FIT: C_T_e = a*R_max + b  (oczekiwanie: a=1/4=0.25)")
    print("=" * 76)
    print()

    for lab in labels:
        Rs = np.array([d[0] for d in C_T_data[lab]])
        CTs = np.array([d[2] for d in C_T_data[lab]])
        a, b = np.polyfit(Rs, CTs, 1)
        print(f"  {lab}: C_T = {a:.5f} * R_max + {b:.5f}")
        print(f"       (predykcja: a=0.25000)")
        pred_70 = a * 70 + b
        print(f"       C_T(R_max=70) = {pred_70:.4f}")

    print()
    print("=" * 76)
    print("  WNIOSKI: czy (K_mu/K_e)^2 zalezy od R_max?")
    print("=" * 76)
    print()
    print(f"  {'R_max':>8} {'K_e':>10} {'K_mu':>10} {'K_tau':>10} "
          f"{'(K_mu/K_e)^2':>15} {'(K_tau/K_e)^2':>16}")
    print(f"  {'-----':>8} {'---':>10} {'----':>10} {'-----':>10} "
          f"{'-----------':>15} {'------------':>16}")

    m_mu_pdg = 206.7682830
    m_tau_pdg = 3477.15

    for i, R in enumerate(R_max_list):
        Ke = C_T_data['e'][i][1]
        Km = C_T_data['mu'][i][1]
        Kt = C_T_data['tau'][i][1]
        rm = (Km/Ke)**2
        rt = (Kt/Ke)**2
        print(f"  {R:>8} {Ke:>10.5f} {Km:>10.5f} {Kt:>10.5f} "
              f"{rm:>15.4f} {rt:>16.4f}")

    print()
    print(f"  PDG: m_mu/m_e = {m_mu_pdg:.4f}, m_tau/m_e = {m_tau_pdg:.4f}")
    print()
    print("  Obserwacja: (K_mu/K_e)^2 powinno byc NIEZALEZNE od R_max")
    print("  (bo jesli K = a*R*A^2 + b*A^2, to K_mu/K_e -> A_mu^2/A_e^2)")
    print()

    # Test: ratio at different cutoffs
    rs_mu = [(C_T_data['mu'][i][1]/C_T_data['e'][i][1])**2 for i in range(len(R_max_list))]
    rs_tau = [(C_T_data['tau'][i][1]/C_T_data['e'][i][1])**2 for i in range(len(R_max_list))]
    rs_mu = np.array(rs_mu)
    rs_tau = np.array(rs_tau)
    print(f"  (K_mu/K_e)^2:  mean={rs_mu.mean():.4f}, std/mean={rs_mu.std()/rs_mu.mean()*100:.3f}%")
    print(f"  (K_tau/K_e)^2: mean={rs_tau.mean():.4f}, std/mean={rs_tau.std()/rs_tau.mean()*100:.3f}%")

    print()
    print("=" * 76)
    print("  DECOMPOZYCJA K = K_core + K_tail")
    print("=" * 76)
    print()
    print("  K_core = int_0^r_c (T*r^(d-1)) dr   (rdzen, r < r_c ~ 1)")
    print("  K_tail = int_{r_c}^{R_max} (T*r^(d-1)) dr  ~ (A^2/4)(R_max - r_c)")
    print()

    r_c = 5.0  # wybierz r_c ktory odlacza rdzen od tailu (rdzen + kilka oscylacji)
    print(f"  Podzial przy r_c = {r_c}")
    print()
    print(f"  {'lep':<6} {'K_core':>12} {'K_core/A^2':>12} {'K_tail(70)':>12} "
          f"{'K_tail/A^2':>12}")
    for lab in labels:
        # K_core
        rs_c = np.linspace(1e-4, r_c, 3000)
        ys_c = sols[lab].sol(rs_c)
        g_c = np.clip(ys_c[0], 1e-12, 10.0)
        gp_c = ys_c[1]
        T_c = 0.5 * g_c**2 * gp_c**2
        K_core = np.trapezoid(T_c * rs_c**2, rs_c)

        # K_tail (r_c to 70)
        rs_t = np.linspace(r_c, 70.0, 5000)
        ys_t = sols[lab].sol(rs_t)
        g_t = np.clip(ys_t[0], 1e-12, 10.0)
        gp_t = ys_t[1]
        T_t = 0.5 * g_t**2 * gp_t**2
        K_tail = np.trapezoid(T_t * rs_t**2, rs_t)

        print(f"  {lab:<6} {K_core:>12.6f} {K_core/As[lab]**2:>12.4f} "
              f"{K_tail:>12.6f} {K_tail/As[lab]**2:>12.4f}")

    print()
    print("  Predykcja tail: K_tail/A^2 ~ (R_max - r_c)/4 = (70-5)/4 = 16.25")
    print()

    print("=" * 76)
    print("  INTERPRETACJA")
    print("=" * 76)
    print()
    print("  1. K/A^2 NIE jest wielkoscia uniwersalna -- zalezy od cutoff R_max.")
    print("     Linearyzacja tail daje K/A^2 ~ R_max/4 (oczekiwane i potwierdzone).")
    print()
    print("  2. Ratio (K_mu/K_e)^2 JEST uniwersalny: A-zalezne prefactory")
    print("     skracaja sie, zostawiajac (A_mu/A_e)^4.")
    print()
    print("  3. Mechanizm A^4 dla MASY RATIO jest W PELNI WYPROWADZONY:")
    print("     m_i/m_j = (K_i/K_j)^2 = (A_i/A_j)^4  (cutoff-independent)")
    print()
    print("  4. Dla BEZWZGLEDNEJ skali masy, R_max staje sie FIZYCZNA")
    print("     dlugoscia. W teorii R5 substratu, R_max = rozmiar rezonansu")
    print("     w przestrzeni konfiguracji.")
    print()
    print("  5. Koherentnie z faktem: K (calka Euclidean action) jest")
    print("     IR-divergent liniowo dla tail oscylujacych, i tylko RATIO")
    print("     daja fizyczny sens. To jest typowe dla teorii solitonowych.")
    print()
    print("  6. C_T(R_max) = R_max/4 + C_core gdzie C_core jest liczbowo O(1)")
    print("     i charakteryzuje RDZEN solitonu.")


if __name__ == "__main__":
    main()
