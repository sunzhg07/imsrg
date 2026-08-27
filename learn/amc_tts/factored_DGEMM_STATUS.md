# Tensor ethS vs scalar Factorized — DGEMM status

**Agent recall:** full playbook in [LESSONS.md](LESSONS.md) §Playbook (read before any new diagram).

Product ranks: \(T\times T\to S\), \(T\times S\to T\), \(S\times S\to S\). Never \(T\times T\to T\).

## Verdict

**Closer.** Production paths match TTS for all Γ diagrams (III/IV). Factorized mid-J / Path B kept in-tree but **disabled** (`kUseFactorized* = false`) until calibrated — they already disagree with AMC TTS at **λ=0** (not a mid-J-only bug).

**Oracle note:** TTS / Path A = truth. Factorized ports from `FactorizedDoubleCommutator.cc` are a different partition; do not claim Fact≡TTS until elementwise calibrated.

| Diagram | χ rank | Scalar Factorized | ethS (λ=0 default) | ethS (λ≠0) |
|---|---|---|---|---|
| \(f^{\mathrm{I,II,III}}\) | — | DGEMM / tables | mostly DGEMM | AMC/tables |
| \(\Gamma^{\mathrm{I,II}}\) | 1b | DGEMM | DGEMM | DGEMM / λ=0 only for II |
| \(\Gamma^{\mathrm{III}_a}\) | scalar χ^η | Pandya→inv→ladder DGEMM | **Factorized IIa/IIc DGEMM** | **mid-J Factorized** (may ≠ TTS) |
| \(\Gamma^{\mathrm{III}_b}\) | scalar χ^η | Pandya→RC→DGEMM | **Factorized DGEMM** | **mid-J Factorized** (may ≠ TTS) |
| \(\Gamma^{\mathrm{III}_c}\) | scalar χ^θ | IIe Pandya→DGEMM | **Factorized IIe/IIf DGEMM** | **mid-J Factorized** (may ≠ TTS) |
| \(\Gamma^{\mathrm{IV}_a}\) | tensor χ | CHI_VI DGEMM | Factorized DGEMM | **Path A** |
| \(\Gamma^{\mathrm{IV}_b}\) | tensor χ | CHI_V DGEMM | Factorized DGEMM | **Path A** |
| \(\Gamma^{\mathrm{IV}_c}\) | tensor χ | CHI_VII DGEMM | Factorized DGEMM | **Path A** |

## Flags

| Piece | default | slow / TTS |
|---|---|---|
| GIIIa | Path B (Pandya/inv χ^η → Chi_AS×Γ), all λ | none (no `tts_GIIIa` / `*_slow`) |
| GIIIb | Factorized (all λ; mid-J at λ≠0) | `SetUse_TypeGIIIb_slow(True)` → TTS |
| GIIIc | Factorized (all λ; mid-J χ^θ at λ≠0) | `SetUse_TypeGIIIc_slow(True)` → Path A |
| GIVa | Factorized λ=0 | `*_slow`→TTS; λ≠0 → Path A |
| GIVb | Factorized λ=0 | `*_slow`→TTS; λ≠0 → Path A |
| GIVc | Factorized λ=0 | `*_slow`→TTS; λ≠0 → Path A |

## Bench (emax=1 He4, rough)

| Piece | λ=0 Factorized | λ=0 PathA/TTS | λ=2 |
|---|---|---|---|
| GIIIa | ~1 ms | PathA ~0.3 s / TTS ~6 s | mid-J Fac (≠ TTS OK) |
| GIIIb | ~1 ms | TTS (via slow) | mid-J Fac (≠ TTS OK) |
| GIIIc | ~1 ms | PathA ~20 ms / TTS ~0.8 s | mid-J Fac (≠ TTS OK) |
| GIVa | ~0 ms | TTS (via slow) | PathA ~1 s / TTS ~7 s |
| GIVb | ~0 ms | TTS (via slow) | PathA ~1.3 s / TTS ~74 s |
| GIVc | ~0 ms | TTS (via slow) | PathA ~0.08 s / TTS ~17 s |

## Gaps to close (priority)

1. ~~**III λ≠0 mid-J** (\(T\times T\to S\))~~ — done for III_a/b/c
2. **IV Path B (optional):** rectangular CC DGEMM — Path A already fast enough for emax=1 smoke
3. Optional: calibrate Factorized vs TTS if TTS is production truth
4. Optional: shared scalar/tensor Pandya helpers (dedupe III_a/b mid-J χ^η)