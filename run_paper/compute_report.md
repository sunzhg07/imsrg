# Compute report — Mg24 E1 and O14→N14 GT (2026-11-11)

Interaction for everything below: EM1.8/2.0 (NN `TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz`
+ NO2B 3N `NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin`), ħω = 16 MeV.
Solver as in `he8_mr_eom.py`: HF → Magnus, Hunter–Gatherer (‖Ω‖ ≤ 0.1), `atan` (s_max 50) then
`shell-model-atan` (s_max 100). Hs and operators are repackaged normal-ordered w.r.t. the inert core.

## 1. Mg24 — VS-IMSRG (ref Mg24, valence sd-shell, core O16), emax = 4

Script: `mg24_vs_imsrg.py` (`EMAX=4`).

| quantity | MeV |
|---|---|
| E_HF (Mg24 ensemble ref) | −115.465 |
| after core decoupling (`atan`) | −136.150 |
| after valence decoupling (`shell-model-atan`) | −155.981 |
| E_core (O16, core-NO) | −113.077 |

16 "BCH_Transform not converged after 12 nested commutators" warnings early in the valence stage
(‖η₂‖ ≈ 5.6 at start); the flow converges (‖η₂‖ → 1e-3, E0 stable to 1e-5). Same warning appears in
`run2/ca40_fp_vs_2b3b.log`.

Outputs (`hs_cache/`): `Hs_Mg24_sdshell_EM1.8_2.0_hw16_emax4.{dat,snt}`,
`E1c_…emax4.dat` (recoil-corrected E1, e_p = N/A, e_n = −Z/A, evolved with the same Ω),
`E1p_…emax4.dat` (bare proton E1, evolved).
Evolved E1c: ‖1b‖ 23.10 → 22.32, ‖2b‖ = 27.49.

## 2. Mg24 — KSHELL sd-shell (4p + 4n) and reference densities

Script: `mg24_kshell_rdm.py`. Work dir `kshell_mg24/e4/`; refs `run2/mg24/mg24_{gs,exc0,2p_gs}_e4.ref`,
summary `run2/mg24/mg24_kshell_e4.json`.

| state | E_SM | E_SM + E_core | E_x | exp. E_x |
|---|---|---|---|---|
| 0⁺₁ | −67.316 | **−180.393** | 0 | 0 (B.E. −198.26) |
| 2⁺₁ | −65.206 | −178.282 | 2.11 | 1.37 |
| 0⁺₂ | −57.694 | −170.771 | 9.62 | 6.43 |
| 0⁺₃ | −55.871 | −168.947 | 11.45 | — |

⟨H⟩_ρ from each density reproduces E_SM + E_core to < 3e-6 MeV, both via `LoadKshellRdm` and via
reading the `.ref` back with `EOM(Hs, ref, …)` (agreement 2e-9 after the `WriteTdm` fix, §5).
emax = 4 is far from converged (g.s. ≈ 18 MeV underbound, spectrum too spread).

## 3. Mg24 — 0⁺ → 1⁻ MR-EOM and E1 strength

Script: `mg24_e1_eom.py`. Reference 0⁺₁ (`run2/mg24/mg24_gs_e4.ref`), tensor EOM J^π = 1⁻
(dimension 16 414). For each Ritz state χ_ν:

  T_ν = `NormMultiref_tensor(E1, χ_ν)` = ½⟨[E1, χ_ν − χ_ν†]⟩_ρ (TT→S commutators, then ⟨·⟩_ρ),
  N_ν = `ComputeNorm(χ_ν, χ_ν)`,
  **B(E1; 0⁺ → 1⁻_ν) = −λ̂² T_ν² / N_ν** (λ̂² = 3).

Normalization check (`TEST=1`, bare E1, single particle–hole χ = (a†ĩ)^{1⁻}): T = ⟨a‖E1‖i⟩·N exactly in all
cases; for a full core hole → empty orbit (0p₃/₂ → 1d₅/₂, physical reduced norm 1) the code gives
N = −1/3, i.e. code N, T = −(1/λ̂²) × physical.

| pair | N (code) | T (code) | ⟨a‖E1‖i⟩ |
|---|---|---|---|
| 0p₃/₂ core → 0d₅/₂ val | −4.924e-2 | −9.487e-2 | 1.9268 |
| 0p₁/₂ core → 1s₁/₂ val | −2.142e-1 | −1.376e-1 | 0.6423 |
| 0d₅/₂ val → 0f₇/₂ q | −2.841e-1 | −7.741e-1 | 2.7250 |
| 1s₁/₂ val → 1p₃/₂ q | −1.191e-1 | −1.710e-1 | 1.4362 |

**Open issue — tensor metric.** With the KSHELL occupations n(0d₅/₂) = 2.95/6 = 0.49 and
n(1s₁/₂) = 0.41/2 = 0.21, the physical reduced norms are 1 − n_a (core → a) and n_i (i → q). The code
values correspond to n = 0.852 and 0.357, i.e. √3 × the true occupations, consistently from both pair
types. So in `ComputeNorm` / `NormMultiref_tensor` the density-dependent terms scale as −1/λ̂ while the
core term scales as −1/λ̂² — one of them is off by λ̂ (invisible for λ = 0). The E1 starting vector
(projected E1 excitation part) has N = −1.93 after `ProjectOprator`, and the He8 tensor runs stop on
`indefinite_metric`; both are consistent with mis-weighted metric terms.
Next: m-scheme check of the tensor norm for single ph / qv configurations against the same density,
to find which term carries the extra λ̂, before changing `EOM.cc`.

Status of the E1 run: started from a random projected vector (as `EOM::RunMR`), MAX_ITER = 60,
STATE_WANT = 30 — **cancelled locally after setup (no results); rerun on riker**; results go to
`eom_result/Mg24_EM1.8_2.0_hw16_emax4_e1_1m.json`. B(E1) values are provisional until the metric
issue is resolved.

Rerun: `MAX_ITER=60 STATE_WANT=30 python3 mg24_e1_eom.py` (`SEED=e1` to start from the E1 vector,
`TEST=1` for the normalization check).

## 4. O14 (0⁺) → N14 (1⁺) Gamow–Teller, p-shell, ref N14, emax = 3

Script: `o14_n14_gt.py`. One VS-IMSRG decoupling (ref N14, valence p-shell, core He4) for both nuclei:
E_HF = −67.549, E0 = −87.759, E_core(He4) = −8.430 MeV. Evolved GT ⟨f‖στ_±‖i⟩ (no g_A):
‖1b‖ 38.52 → 38.34, ‖2b‖ = 22.14 (`hs_cache/GT_N14_pshell_EM1.8_2.0_hw16_emax3.dat`).
KSHELL: O14 = 6p + 4n, N14 = 5p + 5n valence; work dir `kshell_o14n14/e3/`.

KSHELL bare B(GT) (`transit` O14 → N14):

| N14 state | E_x (MeV) | ⟨1⁺‖στ‖0⁺⟩ | B(GT) |
|---|---|---|---|
| 1⁺₁ | 0 | −0.887 | 0.786 |
| 1⁺₂ | 4.95 | −1.982 | 3.926 |
| 1⁺₃ | 14.75 | −1.107 | 1.226 |
| 1⁺₄ | 21.89 | −0.248 | 0.062 |
| sum | | | **6.000** = 3(Z−N), Ikeda |

Exp.: g.s. B(GT) ≈ 1e-3 (anomalously suppressed), 1⁺ at 3.95 MeV ≈ 2.8.

Not done: B(GT) with the evolved operator — needs KSHELL charge-changing OBTD/TBTD, which this
`transit` setup does not print.

EOM note: ⟨N14 1⁺_ν|στ|O14⟩ = ⟨O14|[Q_ν†, στ]|O14⟩ is exactly zero for the bare στ, because the EOM
space has no valence→valence one-body block (⟨q|v⟩, (v∪q)←core, …) and στ only connects the same
n,l (p-shell ↔ p-shell, or 0s core ↔ core). The low N14 1⁺ states are in the shell-model space. The
EOM with the evolved GT gives only the beyond-valence (high-lying) strength — planned after the
tensor-metric issue, using ref O14 0⁺ and a charge-changing 1⁺ EOM.

## 5. Code fixes made (in this commit)

- `src/ReadWrite.cc` `ReadOperator` / `ReadOperatorHuman`: the lower triangle of a tensor operator's
  one-body part is now filled with the phase (−1)^{j_i−j_j} (was missing; scalars unaffected).
  E1 write/read round trip 41 → 2e-5.
- `src/EOM.cc` `WriteTdm`: TBTD J columns written undoubled, matching `ReadOsmRdm` and the external
  OSM `.ref` files (was `J*2`). Every `WriteTdm` file written before this is unreadable or misread:
  `run2/he8/he8_*_e{3,4,5}.ref`, `run2/be8/be8_0p_gs.ref`, `run2/be8/be8_2p_gs.ref` — regenerate from
  their transit logs. No EOM result used them (He8 runs load transit logs via `LoadKshellRdm`).

## 6. To do on riker

1. Build at this commit; rerun `mg24_vs_imsrg.py` (also writes E1c/E1p), `mg24_kshell_rdm.py`.
2. m-scheme check of the tensor norm (§3 open issue), fix, then `mg24_e1_eom.py`.
3. Mg24 at emax = 6 (`EMAX=6` for all three scripts).
4. Regenerate He8/Be8 `WriteTdm` refs.
5. GT with the evolved operator (KSHELL charge-changing densities); EOM beyond-valence GT.
