# Neutrina i oscylacje w gęstym ośrodku (efekt MSW z TGP)

**Data startu:** 2026-04-20
**Status:** plan / open
**Kategoria:** TGP applications → leptony → oscylacje neutrin

## 1. Problem fizyczny

**Efekt Mikheyeva–Smirnova–Wolfensteina (MSW)** (Wolfenstein 1978, Mikheyev-Smirnov 1985):
neutrino oddziałujące ze zwykłą materią nabywa efektywną masę
$$A = 2\sqrt{2}\,G_F N_e E_\nu$$
gdzie $N_e$ = gęstość elektronów. Hamiltonian oscylacji (2-flavor):
$$H = \begin{pmatrix}-\Delta_m^2/4E\cos 2\theta + A/2 & \Delta_m^2/4E\sin 2\theta \\
\Delta_m^2/4E\sin 2\theta & \Delta_m^2/4E\cos 2\theta - A/2\end{pmatrix}$$

MSW wyjaśnia rezonansową konwersję $\nu_e\to\nu_\mu$ w Słońcu (rozwiązanie
problemu neutrin słonecznych, SNO 2001) i wpływa na neutrina z supernowych.

**Otwarte problemy:**
- masy bezwzględne neutrin (hierarchia, NH vs IH — T2K/NOνA wciąż niejednoznaczne),
- faza CP Diraca $\delta_\text{CP}$ — dopiero się mierzy,
- "non-standard interactions" (NSI): uogólnione sprzężenia neutrina-materia,
  które mogą imitować MSW + modyfikować $\theta_{13}, \theta_{23}$.

## 2. Dlaczego TGP

W TGP **przestrzeń jest generowana przez masę**: gęsta materia → lokalnie
większe $\Phi$, silniejszy gradient, zmodyfikowana metryka efektywna
$g_{ij} = e^{+2U/c_0^2}\delta_{ij}$.

**Hipoteza kluczowa:** macierz mieszania neutrin PMNS nie jest fundamentalna,
lecz jest **efektywna** — w próżni Φ=Φ₀ daje kąty wakuumowe, a w gęstej
materii "obraca się" razem z gradientem $\Phi$:
$$U_{\text{PMNS}}^{\text{eff}}(x) = U_0 \cdot R\big(\theta_\Phi(\nabla\Phi)\big).$$

W tym sensie efekt MSW **nie wymaga oddziaływania W-boson z elektronami**,
lecz jest czysto geometryczną konsekwencją zmiany substratu — pozornie
bliskie "gravitational MSW" (Grossman-Nir 1997) ale wyprowadzone z pierwszej
zasady TGP.

Powiązanie z Koide i masami leptonów (istniejący sektor `research/particle_sector_closure`):
skoro $m_e, m_\mu, m_\tau$ wychodzą ze struktury substratu, masy neutrin
mogą być $m_i = m_i^{(0)}\cdot f(\Phi/\PhiZero)$.

## 3. Cele badawcze

### N1 — MSW z geodezyjnego propagatora

Wyprowadzić Hamiltonian oscylacji z propagatora na krzywej geodezyjnej
substratu:
$$i\partial_\tau\,\psi = \big(m_\text{eff}^2(\Phi)/2E\big)\,\psi$$
gdzie $m_\text{eff}^2(\Phi)$ zawiera klasyczny Wolfensteina $2\sqrt2 G_F N_e E$
jako człon liniowy w materii, ale TGP przewiduje **wyższe człony**
$\propto (\nabla\Phi)^2/c_0^4$ — testowalne w ekstremalnych gęstościach.

### N2 — Słoneczne neutrina z TGP

Używając rzeczywistego profilu $\rho_\text{Sun}(r)$ i $\Phi(r)$ z ansatzu
`research/metric_ansatz`, obliczyć $P(\nu_e\to\nu_e)$ i porównać z danymi
SNO, Borexino, Super-K. Sprawdzić, czy 3-flavor fit daje ten sam
$\theta_{12}, \theta_{13}, \Delta m_{21}^2$ co Standard Model w 1σ.

### N3 — Supernowa SN1987A + przyszłe SN

Dla gęstości $\rho\sim10^{14}$ g/cm³ w proto-neutron star, $\Phi$ jest bliskie
silnej deformacji. Przewidywać modyfikację **neutrino sphere** i widma
neutrin. Dopasowanie do Kamiokande/IMB 1987 + gotowość na następną
galaktyczną SN (JUNO, DUNE, Hyper-K).

### N4 — Predykcja masy bezwzględnej

Jeśli $m_\nu$ jest generowana przez minimalną perturbację $\Phi$ od
elektroosłabego próżniowego $\langle H\rangle$, TGP może dać
$\sum m_\nu$ z $\beta, \gamma, \PhiZero$ (a nie fittować).
Cel: zawęzić do $\sum m_\nu = (40\pm 10)$ meV i porównać z KATRIN +
kosmologią ($<0.12$ eV, Planck 2018).

## 4. Plan numeryczny

- **ps01_MSW_vacuum_baseline.py** — 3-flavor propagacja w próżni.
- **ps02_MSW_solar_standard.py** — klasyczny MSW w Słońcu, reprodukcja
  PDG best-fit.
- **ps03_MSW_TGP_substrate.py** — dodatkowy człon $(\nabla\Phi)^2/c_0^4$
  i fit do danych SNO/Borexino.
- **ps04_supernova_neutrino_TGP.py** — propagacja przez PNS z TGP-metryką.
- **ps05_absolute_mass_from_substrate.py** — predykcja $\sum m_\nu$ od
  nonlinearity $\Phi^3$ w potencjale.

## 5. Literatura startowa

- Wolfenstein 1978, Mikheyev-Smirnov 1985 — oryginały
- PDG 2024 review *Neutrino masses and mixing*
- SNO: Ahmad et al., PRL 87, 071301 (2001) — flavor conversion confirmation
- KATRIN: Aker et al., *Direct neutrino-mass measurement* Nature Phys. 18 (2022), m < 0.8 eV
- Grossman-Nir 1997 — gravitational contribution to MSW
- Gonzalez-Garcia-Maltoni, Phys. Rept. 460, 1 (2008) — NSI review

## 6. Relacje z innymi sektorami TGP

- **particle_sector_closure + brannen_sqrt2**: masy leptonów $m_e, m_\mu, m_\tau$
  z substratu. Ta sama struktura może generować $m_{\nu_i}$.
- **cabibbo_correction**: kąty mieszania kwarków z TGP → analog dla PMNS.
- **metric_ansatz**: profile $\Phi$ wewnątrz gwiazd/SN.
- **cosmo_tensions**: sum $m_\nu$ wchodzi w CMB+LSS (S₈).

## 7. Falsyfikowalność

- Jeśli TGP daje konkretne $c_2\ne 0$ w członie $(\nabla\Phi)^2$, to widmo
  dzień/noc Super-K dla neutrin słonecznych powinno odchylać od czystego
  MSW — mierzalne do <1% na obecnych statystykach.
- Predykcja $\sum m_\nu \in (30,50)$ meV byłaby testowalna w DESI + CMB-S4
  w ciągu 5 lat.

## 8. Otwarte pytania

1. Jak zdefiniować "efektywną masę" neutrina w TGP-substratie bez zakładania,
   że neutrino jest Diracowskie vs Majoranowskie?
2. Czy $\delta_\text{CP}$ wynika z chiral Z_2 symmetry substratu?
3. Czy "neutrino-TGP coupling" narusza lepton number / baryon asymmetry?

## 9. Link do rdzenia TGP

Core paper [[tgp_core.pdf]] Theorem *Matter generates space*
(A-IV) i § *Emergent metric* — podstawa dla geometrycznego wyprowadzenia
oscylacji.

Status w TGP status-tree: nowa gałąź *leptons → neutrinos*, uzupełniająca
`particle_sector_closure` (masy naładowanych leptonów).
