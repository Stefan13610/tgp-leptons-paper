# Theory of Generated Space — Particle-sector closure (charged leptons)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19706861.svg)](https://doi.org/10.5281/zenodo.19706861)

This repository accompanies the preprint:

> **Particle-sector closure: a TGP-structured derivation of charged-lepton masses, the Koide relation, and the Cabibbo angle**
> M. Serafin, 2026.
> Zenodo DOI [10.5281/zenodo.19706861](https://doi.org/10.5281/zenodo.19706861).

It is a companion to the TGP core paper:

> **Theory of Generated Space: A minimal core of axioms, substrate, and effective field** —
> Zenodo DOI [10.5281/zenodo.19670324](https://doi.org/10.5281/zenodo.19670324),
> repository [Stefan13610/tgp-core-paper](https://github.com/Stefan13610/tgp-core-paper).

Starting from the substrate-field framework fixed by the core paper, we
derive, from the same soliton ODE, (i) the charged-lepton mass law
$m_n = c_M\,A_{\text{tail},n}^{4}$, (ii) the Koide relation $K=2/3$ as a
theorem of the Brannen ratio $B=b/a=\sqrt{2}$, (iii) the integer $N=3$
generation count as a metric-singularity barrier, and (iv) the Cabibbo
angle $\lambda_C$ through a $\mathbb{Z}_3$ self-energy subtraction on
$\mathrm{GL}(3,\mathbb{F}_2)$. The residual analytic primitive
$P_{\cos}$ in $\alpha_3$ is declared and bounded.

**Headline numerical results** (verified, see `research/`):

- Mass ratio $m_\mu/m_e = (A_\mu/A_e)^4 = 206.74$ vs PDG $206.768$ — **0.013 %**.
- Koide $K = 2/3$ derived: $g_0^\tau$ from mass matches $g_0^\tau$ from
  Koide to within **0.036 %** (1.73027 vs 1.72965).
- $B = \sqrt{2}$ holds to $|B_{\text{num}} - \sqrt 2| < 10^{-6}$ over
  $N=50\,000$ scan of PDG-consistent triples.
- $\alpha_2 = (1 - \ln 3 / 4)/2$ — analytic (Frullani).
- $\alpha_3 = \pi^2/128 + P_{\cos}$ — structural decomposition;
  $P_{\cos}=0.012615939290\ldots$ at 40-digit precision, no relation to
  any 24-element polylog/Clausen basis at $\max\mathrm{coef}=10^{14}$,
  promoted to TGP-named constant `TGP_Pcos`.
- Cabibbo angle $\lambda_C = (0.6847/3)\cdot 165/167 = 0.22550$ vs PDG
  $0.22500$ — **0.75 $\sigma$** (from 4.8 $\sigma$ at zero order).
- Quark Koide shown RG-invariant under 2-loop QCD $\gamma_m$
  ($K_u = 0.8746$, $K_d = 0.7398$) — Koide is charged-lepton-specific.
- Neutrino Koide $K_\nu \le 0.553$ (NO) excludes the Brannen $\sqrt{2}$
  ansatz for neutrinos by oscillation $\Delta m^2$ data.

**Open problems explicitly flagged**: analytic identification of
$P_{\cos}$ in a basis larger than 30 polylog/Clausen constants.

## Repository contents (planned layout)

```
paper/
  tgp_leptons.tex                — LaTeX source
  tgp_leptons.pdf                — compiled preprint

research/                        — numerical support to be cited in the paper
  particle_sector_closure/       — ps1..ps4: α₃, g₀^τ, quark Koide, neutrino Koide
  brannen_sqrt2/                 — B=√2, Koide, c₁ = 1−ln3/4, α_n perturbative
  why_n3/                        — N=3 from g₀_crit barrier + A⁴ mass law
  mass_scaling_k4/               — m = c · K² virial mechanism, k=4 uniqueness
  cabibbo_correction/            — GL(3,𝔽₂) Z₃ self-energy subtraction
  neutrino_msw/                  — (PLAN) Majorana/seesaw mass-ratio structure
```

## Extraction plan (what to import from the working workshop)

Source tree for the snapshot is
[Stefan13610/TGP](https://github.com/Stefan13610/TGP)
(`TGP_v1/research/`). The following files are the minimal closure set
and will be copied verbatim into `research/` on release:

### From `particle_sector_closure/` (the master consolidator, P4 closure 3/4 full + 1/4 structural)

- `README.md` (the P4 closure document)
- `ps1_alpha3_hunting.py` + `ps1_results.txt` — α₃ at dps=40, 24-element PSLQ basis to max_coef=10¹⁴
- `ps2_g0_tau_ODE_scan.py` + `ps2_results.txt` — g₀^τ from A⁴ vs Koide bisection (0.036 % agreement)
- `ps3_quark_koide_qcd.py` + `ps3_results.txt` — 2-loop QCD running, K_u, K_d RG-invariant
- `ps4_neutrino_majorana.py` + `ps4_results.txt` — neutrino Koide + Brannen ansatz falsification

### From `brannen_sqrt2/` (analytic Brannen / Koide backbone)

Required scripts (not the entire 80-script exploratory tree):
- `r6_c1_closed_form_test.py`, `r6_c1_perturbative.py` — c₁ = 1−ln3/4
- `r6_c11_alpha3_ultra_precision.py` — α₃ decomposition at ultra-HP
- `r6_c15_alpha3_pslq_identify.py`, `r6_c17_pslq_comprehensive.py`, `r6_c18_pslq_extended.py` — PSLQ scans
- `r6_koide_from_ode.py`, `r6_koide_variational.py` — Koide from soliton ODE
- `r6_eta_koide_attack.py` — eta-level consistency with Koide
- `r6_tau_constraint.py` — τ constraint on g₀
- `README.md`

Bundle the ~70 exploratory `r6_c2_*` / `r6_c2[0-9]_*` scripts as `research/brannen_sqrt2/exploratory/` without individual citation.

### From `why_n3/`

- `README.md`
- `r3_n3_from_barrier.py`, `r3_g0crit_analytical.py`, `r3_self_space_stability.py` — g₀_crit = 2.206 and the barrier argument
- `r3_mass_A4_derivation.py`, `r3_mass_function.py`, `r3_virial_mechanism.py` — m = A⁴ mechanism
- `r3_koide_derivation.py`, `r3_koide_pi_over_k.py` — Koide from geometry
- `r3_conservation_universal.py`, `r3_sum_conservation.py` — conservation laws
- `r3_alpha_scan.py`, `r3_CT_analytical.py` — α-dependence and topological core

### From `mass_scaling_k4/`

- `README.md`
- `r5_virial_mass_derivation.py`, `r5_mass_ratio_verification.py` — k=4 proof
- `r5_k_squared_mechanism.py` — m = c · K² non-perturbative mechanism (7/7 PASS)
- `r5_e3_cancellation.py` — why perturbative E^(3) does not vanish (context)

### From `cabibbo_correction/`

- `README.md`
- `r1_cabibbo_correction_derivation.py` — Z₃ subtraction formula, 0.75 σ from PDG
- `r1_gl3f2_structure.py` — GL(3,𝔽₂) group structure

### From `neutrino_msw/`

- `PLAN.md` only (the Majorana + MSW programme is a follow-up; neutrino
  Koide falsification in `ps4` is sufficient for this paper)

## Build

```
cd paper
pdflatex tgp_leptons.tex
pdflatex tgp_leptons.tex        # second pass for cross-references
```

Requirements: any modern LaTeX distribution with `amsmath`, `amssymb`,
`amsthm`, `mathtools`, `longtable`, `enumitem`, `hyperref`, `caption`.

## Reproducibility of the numerical support

All Python scripts under `research/` run standalone:

```
python research/particle_sector_closure/ps1_alpha3_hunting.py
python research/particle_sector_closure/ps2_g0_tau_ODE_scan.py
python research/particle_sector_closure/ps3_quark_koide_qcd.py
python research/particle_sector_closure/ps4_neutrino_majorana.py
python research/brannen_sqrt2/r6_c1_perturbative.py
python research/why_n3/r3_n3_from_barrier.py
python research/mass_scaling_k4/r5_k_squared_mechanism.py
python research/cabibbo_correction/r1_cabibbo_correction_derivation.py
```

`ps1` requires `mpmath` at `dps=40`; everything else needs only `numpy`,
`scipy` and `matplotlib`.

## Paper outline (draft)

1. **Introduction** — charged-lepton mass problem, Koide history, TGP
   positioning, core-paper prerequisites.
2. **Substrate soliton and the mass law.**
   Derivation of $m = c_M\,A_{\text{tail}}^4$ from the virial
   $m_{\text{phys}} = c\,K^2$ mechanism; uniqueness of $k=4$ in $d=3$;
   numerical verification $m_\mu/m_e = 206.74$ (0.013 % vs PDG).
3. **$N=3$ generations from metric-singularity barrier.**
   $g_0^{\text{crit}}(d=3,\alpha=1) = 2.206$; $g_0^\tau < g_0^{\text{crit}}$;
   fourth generation forbidden by mass divergence.
4. **Koide $K=2/3$ as a theorem.**
   Brannen $B=\sqrt{2}$ chain; soliton ODE perturbation expansion;
   $c_1 = 1 - \ln 3 / 4$ (Frullani proof); $\alpha_2 = c_1/2$;
   $\alpha_3 = \pi^2/128 + P_{\cos}$ (structural closure); independent
   derivation of $g_0^\tau$ from PDG and from $K=2/3$ matching to
   0.036 %.
5. **The Cabibbo angle.**
   Zero-order $\lambda_C = \Omega_\Lambda / N = 0.22823$ (4.8 σ tension);
   $\mathbb{Z}_3$ self-energy correction on $\mathrm{GL}(3,\mathbb{F}_2)$
   yields $0.22550$ (0.75 σ).
6. **Sector-specificity of Koide.**
   Quark Koide under 2-loop QCD running — RG invariance of
   $K_u, K_d \neq 2/3$; neutrino Koide with full oscillation data —
   Brannen $\sqrt{2}$ ansatz excluded by $\Delta m^2$. Koide is a
   charged-lepton signature of the soliton ODE, not a universal law.
7. **Open problems.**
   Analytic identification of $P_{\cos}$ beyond polylog/Clausen bases;
   Majorana/seesaw mass-ratio structure for neutrinos.

## Relationship to the core paper and the full workshop

This repository is a deliberately narrow slice:

- The **core paper** (axioms, substrate, effective metric, PPN
  parameters, gravitational-wave propagation, 15 proofs + 10 open
  problems) lives at
  [Stefan13610/tgp-core-paper](https://github.com/Stefan13610/tgp-core-paper)
  and is cited here via its Zenodo DOI.
- The **superconductivity closure paper** lives at
  [Stefan13610/tgp-sc-paper](https://github.com/Stefan13610/tgp-sc-paper)
  and is independent of the present work.
- The **full TGP research workshop** — long-form companion manuscript,
  QM-foundations studies, cosmology channels (Hubble, DESI, S8),
  additional sector-closure work, galaxy / mass scaling, UV completion —
  is kept separately at
  [Stefan13610/TGP](https://github.com/Stefan13610/TGP).

That is where development happens. This repository is the stable,
paper-aligned snapshot for the particle-sector closure only.

## Citation

Please cite both DOIs together:

- **This repository (particle-sector closure):** [10.5281/zenodo.19706861](https://doi.org/10.5281/zenodo.19706861)
- **Core paper:** [10.5281/zenodo.19670324](https://doi.org/10.5281/zenodo.19670324)

BibTeX:

```bibtex
@misc{Serafin2026TGPLeptons,
  author       = {Serafin, Mateusz},
  title        = {{Particle-sector closure: a TGP-structured derivation
                   of charged-lepton masses, the Koide relation, and the
                   Cabibbo angle}},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19706861},
  url          = {https://doi.org/10.5281/zenodo.19706861}
}

@misc{Serafin2026TGPCore,
  author       = {Serafin, Mateusz},
  title        = {{Theory of Generated Space: A minimal core of axioms,
                   substrate, and effective field}},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19670324},
  url          = {https://doi.org/10.5281/zenodo.19670324}
}
```

The `.zenodo.json` file contains the machine-readable metadata Zenodo
uses on each release.

## License

The paper text and the accompanying numerical code are released under
[CC BY 4.0](LICENSE).
