# R1: Korekcja Cabibbo — odejmowanie self-energii Z₃

## Problem

Kąt Cabibbo — największe pojedyncze napięcie w TGP:

```
TGP (zerowy rząd):  λ_C = Ω_Λ/N = 0.6847/3 = 0.22823
PDG:                λ_C = 0.22500 ± 0.00067
Napięcie:           4.8σ
```

Kill criterion K2 przeżywa (< 30% off), ale 4.8σ wymaga wyjaśnienia.

## ROZWIĄZANIE (2026-04-14)

### Korekcja: Z₃ self-energy subtraction

**Formuła:**
```
λ_C = (Ω_Λ / N) × (|G| - |Z₃|) / (|G| - 1)
    = (0.6847 / 3) × 165/167
    = 0.22550
```

**Napięcie po korekcji: 0.75σ** (z 4.8σ → 0.75σ)

### Fizyczna motywacja

1. Mieszanie CKM mediowane przez GL(3,𝔽₂) z 168 elementami
2. Elementy Z₃ (3 sztuki: tożsamość + 2 generatory) **zachowują** numer generacji
3. Nie przyczyniają się do mieszania międzygeneracyjnego
4. "Aktywne" kanały mieszania: |G| - |Z₃| = 165 z |G| - 1 = 167 nietrivialnych

**Czynnik korekcyjny:** F = 165/167 = 0.988024

### Efektywna liczba generacji

```
N_eff = N × (|G| - 1) / (|G| - |Z₃|) = 3 × 167/165 = 3.03636
```

To jest czysto grupowo-teoretyczne N_eff (nie kosmologiczne 3.044).

### Kluczowe wyniki

| Wielkość | Przed korektą | Po korekcie | PDG | σ |
|----------|:---:|:---:|:---:|:---:|
| λ_C | 0.22823 | **0.22550** | 0.22500 ± 0.00067 | **0.75** |
| N_eff | 3.000 | **3.036** | — | — |
| Status | 4.8σ KRYTYCZNE | **0.75σ PASS** | — | — |

### Weryfikacja CKM

| Element | TGP (z korektą) | PDG | σ |
|---------|:---:|:---:|:---:|
| |V_us| | 0.22550 | 0.22500 | 0.7 |
| |V_cd| | 0.22550 | 0.22486 | 1.0 |
| |V_cb| | 0.04200 | 0.04182 | 0.2 |

### Test jednoznaczności grupy

Korekcja 165/167 działa TYLKO dla GL(3,𝔽₂):

| Grupa | |G| | F | λ_C | σ |
|-------|:---:|:---:|:---:|:---:|
| S₃ | 6 | 0.600 | 0.137 | 131 |
| A₄ | 12 | 0.818 | 0.187 | 57 |
| A₅ | 60 | 0.966 | 0.221 | 6.7 |
| **GL(3,𝔽₂)** | **168** | **0.988** | **0.226** | **0.75** |
| GL(4,𝔽₂) | 20160 | 0.9999 | 0.228 | 4.8 |

→ **Potwierdza, że GL(3,𝔽₂) jest jedyną grupą smakową dającą zgodność.**

## Struktura GL(3,𝔽₂) (z analizy)

```
|GL(3,𝔽₂)| = 168 = PSL(2,7)
Klasy sprzężoności: 6
  1A: rozmiar 1,   rząd 1  (tożsamość)
  2A: rozmiar 21,  rząd 2  (inwolucje)
  3A: rozmiar 56,  rząd 3  (generatory Z₃)
  4A: rozmiar 42,  rząd 4
  7A: rozmiar 24,  rząd 7
  7B: rozmiar 24,  rząd 7

Podgrupy Z₃: 28 (jedna klasa sprzężoności)
|N_G(Z₃)| = 6,  N(Z₃)/Z₃ = Z₂ = Aut(Z₃)
Podwójne koklasy Z₃\G/Z₃: 20
```

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `r1_gl3f2_structure.py` | Pełna analiza algebraiczna GL(3,𝔽₂) | ✅ DONE |
| `r1_cabibbo_correction_derivation.py` | Derywacja korekcji 165/167 | ✅ DONE |

## Pliki do scalenia z rdzeniem (DO ZROBIENIA)

- `tgp_companion.tex` §F2: nowy paragraf z korekcją 165/167
- `scripts/cabibbo_correction_verify.py`: skrypt weryfikacyjny do CI
- Ewentualnie: nowy ex*.py w `nbody/examples/`

## Referencje rdzenia

- `tgp_companion.tex`, linie 329–342, 791–792
- `tgp_letter.tex`, linia 135
- `nbody/examples/ex247_cabibbo_omega_lambda.py` (walidacja λ_C)
- `nbody/examples/ex274_cabibbo_correction.py` (skan korekcji)

## Status

- [x] Krok 1: Algebraiczna struktura GL(3,𝔽₂) — 6 klas, 28 Z₃, 20 DC
- [x] Krok 2: Korekcja Z₃ self-energy: F = 165/167
- [x] Krok 3: Weryfikacja: 4.8σ → **0.75σ**
- [ ] Krok 4: Scalenie z rdzeniem (tgp_companion.tex)
- [ ] Krok 5: Formalizacja dowodu (Lean 4?)
