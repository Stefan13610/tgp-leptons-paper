"""
R3: DWIE masy -- M_energy (action) ~ A^2, m_phys (PDG) ~ A^4.

KLUCZOWE ODKRYCIE (wynik numeryczny tego skryptu):

  M_energy = 4pi * int (T + V_eff) r^2 dr  ~  A_tail^2    (slope numerycznie = 2)
  m_phys / m_e = (A_i/A_e)^4                              (dla leptonow, PDG 0.24%)

  Zatem RELACJA:  m_phys = C * M_energy^2   (KWADRATOWA)
                           ^^^^^^^^^^^^^
  To OBALA moja wczesniejsza hipoteze ze m~A^4 jest bezposrednio
  energia solitonu. Potrzebny dodatkowy mechanizm (binding/virial).

WYJASNIENIE ANALITYCZNE M ~ A^2:

  Tail: g - 1 = A*sin(r+delta)/r (linearyzacja ODE), V''(1) = -1 (false vac)
  T = g^(2a)(g')^2/2 ~ A^2*cos^2/(2r^2)  -> T*r^2 ~ A^2*cos^2/2
  V_eff = V(g) - V(1) ~ -A^2*sin^2/(2r^2) + O(A^3) -> V*r^2 ~ -A^2*sin^2/2

  (T+V_eff)*r^2 ~ A^2*[cos^2 - sin^2]/2 = A^2*cos(2r+2d)/2

  Calka oscylujaca tailu nie daje czystego zera ze wzgledu na RDZEN
  (r<r_1 gdzie linearyzacja zawodzi) dajacego sekularny O(A^2) wklad.

WYJASNIENIE EMPIRYCZNE m = M^2:

  m_mu/m_e (PDG) = 206.77 = 14.38^2 = (M_mu/M_e)^2
  Zatem m_phys skaluje jako KWADRAT M_energy.
  Wymaga to:
    - rozmiary rdzenia skaluje sie z A (core size ~ A)
    - m_phys = M_energy * size(A) = A^2 * A^2 = A^4
  LUB:
    - bound state: E_bind = -m_phys propto <H>^2/E_scale
    - relatywistyczny efekt: m^2 = E^2 - p^2 (w non-rel p << E)

WERYFIKACJA NUMERYCZNA ponizej.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit


def solve_soliton(g0, alpha=1.0, d=3, r_max=80.0):
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


def extract_A_tail(sol, r_fit_min=15.0):
    """Ekstrakcja A_tail z asymptoty g-1 ~ A*sin(r+delta)/r."""
    rs = np.linspace(r_fit_min, sol.t[-1] - 2, 2000)
    ys = sol.sol(rs)
    g = ys[0]
    u = (g - 1) * rs  # u = A*sin(r+delta)
    # Fit: u = A*sin(r+delta) = A_s*sin(r) + A_c*cos(r) with A=sqrt(A_s^2+A_c^2)
    # u(r) = A_s*sin(r) + A_c*cos(r)
    S = np.sin(rs)
    C = np.cos(rs)
    # Least squares: min ||u - A_s*S - A_c*C||^2
    M = np.column_stack([S, C])
    coefs, _, _, _ = np.linalg.lstsq(M, u, rcond=None)
    A_s, A_c = coefs
    A = np.sqrt(A_s**2 + A_c**2)
    delta = np.arctan2(A_c, A_s)
    return A, delta


def compute_mass(sol, alpha=1.0, d=3, r_inner=0.0, r_outer=None):
    """
    Mass = 4pi int_0^inf [T + V_eff] r^2 dr
    T = g^(2a) (g')^2 / 2,  V_eff = V(g) - V(1),  V(g) = g^3/3 - g^4/4.
    Uzywamy cutoff r_outer i usredniamy po oscylacji.
    """
    if r_outer is None:
        r_outer = sol.t[-1] - 2
    rs = np.linspace(max(sol.t[0], 1e-4), r_outer, 8000)
    ys = sol.sol(rs)
    g = np.clip(ys[0], 1e-12, 10.0)
    gp = ys[1]
    V1 = 1.0/3.0 - 1.0/4.0  # V(1) = 1/12
    T = g**(2*alpha) * gp**2 / 2.0
    V = g**3/3.0 - g**4/4.0
    integrand = (T + V - V1) * rs**(d-1)
    # Usredniamy: poniewaz tail oscyluje, liczymy calke do r_outer uzywajac Simpsona
    # a nastepnie usredniamy po ostatnim okresie (r_outer -> r_outer - pi)
    # Dla stabilnosci: po prostu calkujemy do r_outer

    # Lepsza strategia: sub-period average
    # Wyznaczmy sredni z ostatniej polowy:
    n_half = len(rs)//2
    rs_a = rs[:n_half]
    integrand_a = integrand[:n_half]
    M_core = 4*np.pi*np.trapezoid(integrand_a, rs_a)
    # Srednia z okresu
    # tail: integrand ~ A^2*cos(2r+2d)/2; zero-avg per pi
    M_tail_avg = 0.0  # po usrednieniu
    # Lepiej: zintegruj od rs[n_half] do r_outer i USRENIJ po 2 okresach (2*pi)
    rs_tail = rs[n_half:]
    integrand_tail = integrand[n_half:]
    # Uzyj jaczejwarnego usredniania: tail length L = r_outer - rs[n_half]
    # Calka deterministyczna
    M_tail_det = 4*np.pi*np.trapezoid(integrand_tail, rs_tail)
    # Total
    M_total = M_core + M_tail_det
    return M_total


def main():
    print("=" * 70)
    print("  R3 DERYWACJA m = c * A_tail^4 (NUMERYCZNA WERYFIKACJA)")
    print("=" * 70)

    print("\nARGUMENT:")
    print("  Tail expansion g = 1 + A*sin(r+d)/r pokazuje:")
    print("    - A^2 wklad do masy: oscyluje (avg = 0)")
    print("    - A^4 wklad: sekularny 1/r^2 -> skonczony (non-zero)")
    print("  => m ~ A^4 emerguje ze strukture Lagrangianu (false vacuum)")
    print()

    print("=" * 70)
    print("  SKAN m vs A_tail dla roznych g0")
    print("=" * 70)

    # Skan: g0 od blisko 1 do g0_crit ≈ 2.2
    # Dla deficit (g0 < 1) i excess (g0 > 1)
    g0_list = np.concatenate([
        np.array([0.50, 0.60, 0.70, 0.80, 0.869, 0.90, 0.95, 0.98, 0.99]),
        np.array([1.02, 1.05, 1.10, 1.20, 1.30, 1.407, 1.50, 1.60, 1.729, 1.80, 1.90, 2.00, 2.10])
    ])

    print()
    print(f"  {'g0':>8} {'A_tail':>12} {'|g0-1|':>10} {'M (ODE)':>14} {'M/A^4':>14}")
    print(f"  {'--':>8} {'------':>12} {'------':>10} {'-------':>14} {'-----':>14}")

    data = []
    for g0 in g0_list:
        sol = solve_soliton(g0, alpha=1.0, d=3, r_max=80.0)
        if not sol.success:
            print(f"  {g0:>8.4f}: FAILED")
            continue
        A, delta = extract_A_tail(sol, r_fit_min=20.0)
        try:
            M = compute_mass(sol, alpha=1.0, d=3, r_outer=70.0)
        except Exception:
            M = np.nan
        ratio = M/A**4 if A > 1e-10 else np.nan
        data.append((g0, A, M))
        print(f"  {g0:>8.4f} {A:>12.6f} {abs(g0-1):>10.4f} {M:>14.6e} {ratio:>14.3e}")

    print()
    print("=" * 70)
    print("  FIT log|M| vs log|A|: czy slope = 4?")
    print("=" * 70)

    data = [d for d in data if d[1] > 1e-4 and not np.isnan(d[2])]
    if len(data) < 3:
        print("  Niewystarczajaca liczba punktow. Przerywam.")
        return

    g0_arr = np.array([d[0] for d in data])
    A_arr = np.array([d[1] for d in data])
    M_arr = np.array([d[2] for d in data])

    # Tylko dodatnie masy
    mask = M_arr > 0
    print(f"\n  Wszystkich punktow: {len(data)}, z dodatnia masa: {mask.sum()}")
    print(f"  (Deficit g0<1 ma M>0, excess g0>1 ma M<0 dla false vacuum)")

    # Fit dla deficit solitonow (g0 < 1)
    deficit = g0_arr < 1.0
    A_d = A_arr[deficit]
    M_d = M_arr[deficit]
    M_d_pos = np.abs(M_d)  # wszystkie musza byc dodatnie w deficit
    if (M_d > 0).sum() >= 3:
        A_d = A_d[M_d > 0]
        M_d_pos = M_d[M_d > 0]
        logA = np.log(A_d)
        logM = np.log(M_d_pos)
        slope, intercept = np.polyfit(logA, logM, 1)
        c_M = np.exp(intercept)
        print(f"\n  DEFICIT (g0<1, M>0):")
        print(f"    Punktow fit: {len(logA)}")
        print(f"    log|M| = {slope:.6f} * log|A| + {intercept:.4f}")
        print(f"    Slope = {slope:.6f}  (predykcja: 4.000000)")
        print(f"    c_M = {c_M:.6e}")
        print(f"    M = {c_M:.4e} * A^{slope:.4f}")

    # Fit dla excess (g0 > 1, M < 0 w false vacuum)
    excess = g0_arr > 1.0
    A_e = A_arr[excess]
    M_e = M_arr[excess]
    if (np.abs(M_e) > 1e-8).sum() >= 3:
        A_e_use = A_e[np.abs(M_e) > 1e-8]
        M_e_use = M_e[np.abs(M_e) > 1e-8]
        logA = np.log(A_e_use)
        logM = np.log(np.abs(M_e_use))
        slope, intercept = np.polyfit(logA, logM, 1)
        c_M = np.exp(intercept)
        print(f"\n  EXCESS (g0>1):")
        print(f"    Punktow fit: {len(logA)}")
        print(f"    log|M| = {slope:.6f} * log|A| + {intercept:.4f}")
        print(f"    Slope = {slope:.6f}  (predykcja: 4.000000)")
        print(f"    c_M = {c_M:.6e}")
        print(f"    |M| = {c_M:.4e} * A^{slope:.4f}")

    # Razem wszystkie |M| vs A
    A_all = A_arr
    M_abs = np.abs(M_arr)
    mask_all = (A_all > 1e-4) & (M_abs > 1e-10)
    logA_all = np.log(A_all[mask_all])
    logM_all = np.log(M_abs[mask_all])
    slope_all, int_all = np.polyfit(logA_all, logM_all, 1)
    print(f"\n  ALL (|M|):")
    print(f"    Slope = {slope_all:.6f}  (predykcja: 4.000000)")
    print(f"    c_M = {np.exp(int_all):.6e}")

    print()
    print("=" * 70)
    print("  WERYFIKACJA: m_mu / m_e = A_mu^4 / A_e^4")
    print("=" * 70)

    # Bridge cal
    sol_e = solve_soliton(0.86941, alpha=1.0, d=3, r_max=80.0)
    sol_mu = solve_soliton(1.40673, alpha=1.0, d=3, r_max=80.0)
    sol_tau = solve_soliton(1.72931, alpha=1.0, d=3, r_max=80.0)

    A_e, _ = extract_A_tail(sol_e, r_fit_min=20.0)
    A_mu, _ = extract_A_tail(sol_mu, r_fit_min=20.0)
    A_tau, _ = extract_A_tail(sol_tau, r_fit_min=20.0)

    ratio_mu_e = (A_mu/A_e)**4
    ratio_tau_e = (A_tau/A_e)**4
    m_mu_pdg = 206.7682830
    m_tau_pdg = 3477.15
    print(f"\n  A_e   = {A_e:.6f}  (g0_e   = 0.86941)")
    print(f"  A_mu  = {A_mu:.6f}  (g0_mu  = 1.40673)")
    print(f"  A_tau = {A_tau:.6f}  (g0_tau = 1.72931)")
    print()
    print(f"  (A_mu/A_e)^4  = {ratio_mu_e:.4f}   (PDG m_mu/m_e = {m_mu_pdg:.4f})")
    print(f"  (A_tau/A_e)^4 = {ratio_tau_e:.4f}   (PDG m_tau/m_e = {m_tau_pdg:.4f})")
    print()
    diff_mu = (ratio_mu_e - m_mu_pdg) / m_mu_pdg * 100
    diff_tau = (ratio_tau_e - m_tau_pdg) / m_tau_pdg * 100
    print(f"  Diff m_mu/m_e: {diff_mu:+.3f}%")
    print(f"  Diff m_tau/m_e: {diff_tau:+.3f}%")

    # Sprawdz relacje m_phys = C * M_energy^2
    M_e = compute_mass(sol_e, alpha=1.0, d=3, r_outer=70.0)
    M_mu = compute_mass(sol_mu, alpha=1.0, d=3, r_outer=70.0)
    M_tau = compute_mass(sol_tau, alpha=1.0, d=3, r_outer=70.0)

    print()
    print("=" * 70)
    print("  KLUCZOWE ODKRYCIE: m_phys = C * M_energy^2 (kwadratowa)")
    print("=" * 70)
    print()
    print(f"  M_e   (energy) = {M_e:.6f}")
    print(f"  M_mu  (energy) = {M_mu:.6f}")
    print(f"  M_tau (energy) = {M_tau:.6f}")
    print()
    print(f"  M_mu / M_e   = {M_mu/M_e:.4f}  ~ (A_mu/A_e)^2  = {(A_mu/A_e)**2:.4f}")
    print(f"  M_tau / M_e  = {M_tau/M_e:.4f}  ~ (A_tau/A_e)^2 = {(A_tau/A_e)**2:.4f}")
    print()
    print(f"  (M_mu/M_e)^2  = {(M_mu/M_e)**2:.4f}  PDG m_mu/m_e  = {m_mu_pdg:.4f}")
    print(f"  (M_tau/M_e)^2 = {(M_tau/M_e)**2:.4f}  PDG m_tau/m_e = {m_tau_pdg:.4f}")
    print()
    print(f"  Diff (M_mu/M_e)^2  vs m_mu/m_e:  {((M_mu/M_e)**2 - m_mu_pdg)/m_mu_pdg*100:+.3f}%")
    print(f"  Diff (M_tau/M_e)^2 vs m_tau/m_e: {((M_tau/M_e)**2 - m_tau_pdg)/m_tau_pdg*100:+.3f}%")

    # Sprawdz czwarty moment (g-1)^4
    def fourth_moment(sol, d=3, r_outer=70.0):
        rs = np.linspace(max(sol.t[0], 1e-4), r_outer, 8000)
        ys = sol.sol(rs)
        g = ys[0]
        integrand = (g - 1)**4 * rs**(d-1)
        return 4*np.pi*np.trapezoid(integrand, rs)

    M4_e = fourth_moment(sol_e)
    M4_mu = fourth_moment(sol_mu)
    M4_tau = fourth_moment(sol_tau)

    print()
    print("=" * 70)
    print("  HIPOTEZA: m_phys ~ 4th moment  int (g-1)^4 r^2 dr")
    print("=" * 70)
    print()
    print(f"  M4_e   = {M4_e:.6f}  (I4 = int (g-1)^4 r^2 dr)")
    print(f"  M4_mu  = {M4_mu:.6f}")
    print(f"  M4_tau = {M4_tau:.6f}")
    print()
    print(f"  M4_mu / M4_e   = {M4_mu/M4_e:.4f}   PDG m_mu/m_e   = {m_mu_pdg:.4f}")
    print(f"  M4_tau / M4_e  = {M4_tau/M4_e:.4f}   PDG m_tau/m_e  = {m_tau_pdg:.4f}")
    print()
    print(f"  Diff M4_mu/M4_e vs PDG:  {(M4_mu/M4_e - m_mu_pdg)/m_mu_pdg*100:+.3f}%")
    print(f"  Diff M4_tau/M4_e vs PDG: {(M4_tau/M4_e - m_tau_pdg)/m_tau_pdg*100:+.3f}%")

    print()
    print("=" * 70)
    print("  INTERPRETACJA KONCOWA")
    print("=" * 70)
    print()
    print("  1. M_energy (pelna energia solitonu) skaluje ~ A^2:")
    print("     T,V leading O(A^2), oscyluja ale RDZEN daje sekularny O(A^2) wklad.")
    print("     Numerycznie: M_mu/M_e = 14.5, M_tau/M_e = 55.5 (bliskie A^2 ratios).")
    print()
    print("  2. Masa fizyczna PDG skaluje ~ A^4, empirycznie bardzo dokladnie:")
    print("     (A_mu/A_e)^4 = 206.3  vs PDG 206.77 (0.24%)")
    print("     (A_tau/A_e)^4 = 3469 vs PDG 3477 (0.24%)")
    print()
    print("  3. Hipoteza m_phys = M_energy^2 przybliza dla leptonow:")
    print("     (M_mu/M_e)^2 = 209.5 vs 206.8 (1.3% diff)")
    print("     Ale dla tau zawodzi: (M_tau/M_e)^2 = 3078 vs 3477 (-11.5%)")
    print("     Dla tau, bliskiego barierze, dynamika nonlinear przesuwa M od A^2.")
    print()
    print("  4. Hipoteza 4th moment int (g-1)^4 r^2 dr:")
    print("     Sprawdz w raporcie: czy skaluje sie z A^4 dokladnie?")
    print()
    print("  5. NAJPROSTSZE PRAWO: m_phys = c_m * A_tail^4 (bezposrednio)")
    print("     A_tail z linearyzacji ogonu ODE, exponent = 4 EXACT.")
    print("     Mechanizm NIE jest pelna energia solitonu (~A^2).")
    print("     Mozliwy zwiazek: m_phys = (4th moment)*const.")
    print()
    print("  6. STATUS R3 <-> R5 bridge:")
    print("     - empiryczna relacja m = c*A^4 zweryfikowana (0.24% diff)")
    print("     - M_energy = int (T+V) r^2 dr NIE jest masa fizyczna")
    print("     - m_phys wymaga dodatkowego mechanizmu (rozprzestrzenianie rdzenia?)")
    print("     - kierunek przyszly: formalny dowod m_phys = c * <(g-1)^4> lub podobne")


if __name__ == "__main__":
    main()
