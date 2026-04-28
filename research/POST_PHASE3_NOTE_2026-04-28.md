# Post-Phase-3 audit note (2026-04-28)

This is a short audit note attached to the published TGP-Leptons paper
(Zenodo DOI [10.5281/zenodo.19706861](https://doi.org/10.5281/zenodo.19706861)).

## Verdict: clean after micro-fix — no impact on theorems or conclusions

The 2026-04-28 deep-scan audit (eight check items × four flask papers,
performed against the workshop master `TGP/TGP_v1/` post-Phase-3.R-final
state) found **one numerical micro-drift** in `paper/tgp_leptons.tex`,
and **no other issues**. The fix is purely a numerical refresh of
the muon/electron mass-ratio row; **no theorem, no proof, and no
conclusion is altered**.

## The fix (in-repo, applied 2026-04-28)

| Location | Before | After |
|---|---|---|
| Lines 74–75 (intro) | `(A_μ/A_e)^4 = 206.74`; ratio at `0.013%` | `(A_μ/A_e)^4 = 206.77`; ratio at `0.001%` |
| Line 342 (numerical-results table) | `206.74 & 206.768 & 1.3·10⁻⁴` | `206.77 & 206.768 & 1·10⁻⁵` |
| Line 868 (Conclusion bullet) | `0.013% in the μ/e ratio` | `0.001% in the μ/e ratio` |

The PDG anchor (206.768) is unchanged. The TGP value is now consistent
with the master ledger value carried forward through Phase 1 / 2 / 3
post-closure.

The k=4 selection theorem (Thm. `thm:k4`), the Koide relation (Thm.
`thm:koide-theorem`), the N=3 generation count (Thm. `thm:N3`), and
all other proofs in the paper are **unaffected**.

## What the audit checked

| Check | Result |
|---|---|
| `r_21 = (A_μ/A_e)^4` numerical row | drift identified and fixed |
| `r_31 = (A_τ/A_e)^4` numerical row | clean |
| `k_eff` (k=4 integer selection) | clean |
| `K_koide = 2/3` | clean |
| Cross-references / theorem chain | clean |
| α=2 selection wording (paper convention) | clean (C5-compliant) |
| Path B m_σ²/m_s² = 2 | not claimed by this paper — N/A |
| 4-of-4 UV synthesis | not claimed by this paper — N/A |

## What this note does *not* change

- The deposited paper is **not retracted and not superseded**. It is
  **amended in-repo** (three single-line numerical refreshes at the
  same row of the same table); the underlying theorems and the
  experimental anchor (PDG 206.768) are untouched.
- No new Zenodo deposit is issued for this micro-refresh at this time.
  Whether to issue a v2 deposit later — once a larger batch of
  refreshments accumulates — is left as a downstream decision.

## Cross-reference

For the full Phase 3 cycle results, see:
- `TGP/TGP_v1/research/op-phase3-uv-completion/Phase3_R_final_results.md` (Phase 3 closed 60/60)
- Companion full addendum:
  [`tgp-core-paper/research/POST_PHASE3_ADDENDUM_2026-04-28.md`](../../tgp-core-paper/research/POST_PHASE3_ADDENDUM_2026-04-28.md)
- Master grand total post-Phase-3: **281** (M9 13 + M10 42 + M11 62 + Phase 1 50 + Phase 2 54 + Phase 3 60).
