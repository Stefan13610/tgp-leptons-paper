# Numerical support — particle-sector closure

This directory will hold, on release, the numerical support folders
cited in the paper.  They are extracted from the working repository
[Stefan13610/TGP](https://github.com/Stefan13610/TGP) (`TGP_v1/research/`).

Planned layout:

```
particle_sector_closure/   — P4 closure: ps1..ps4 (α₃, g₀^τ, quark Koide, neutrino Koide)
brannen_sqrt2/             — Brannen B=√2, Koide, c₁, α_n perturbative
why_n3/                    — N=3 from g₀_crit barrier + A⁴
mass_scaling_k4/           — m = c·K² virial mechanism, k=4 uniqueness
cabibbo_correction/        — GL(3,𝔽₂) Z₃ self-energy subtraction
neutrino_msw/              — Majorana / seesaw programme (plan only)
```

See the top-level `../README.md` for the extraction plan (which scripts
are the minimal closure set vs which are exploratory).

Each sub-folder will carry its own `README.md` summarising the
sub-problems, scripts and outcomes, copied verbatim from the source
tree.
