# tensor_pro_final — locked tensor / Path B diagrams

Curated gold pack for the diagrams that have a **m ≡ AMC direct ≡ Path B** chain
(or documented partial). Working notes stay in `factored_*/NOTES.md`; this tree is
the **stable hirearchy** for equations + benches only.

Cross-cutting: [../LESSONS.md](../LESSONS.md) · [../REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md) · [../F_1B_PROCEDURE.md](../F_1B_PROCEDURE.md)

## Hierarchy (every diagram)

```
<diagram>/
  README.md
  01_m_scheme/EQUATION.md     # trusted m-scheme (analyze)
  02_amc_direct/
    EQUATION.md
    input/   → symlink to AMC direct .txt
    output/  → symlink to AMC .tex
  03_pathB/
    EQUATION.md               # Pandya / RME / DGEMM Path B
    input/   → Path B AMC inputs (if any)
    output/  → Path B .tex (if any)
  04_benches/
    README.md                 # how to run the three-way lock
    scripts/ → symlink to run/test_*.py
```

## Diagram index

| Dir | Diagram | Status |
|---|---|---|
| [`fI/`](fI/) | \(f^{\mathrm{I}}\) / χ^α | **DONE** |
| [`fII/`](fII/) | \(f^{\mathrm{II}}\) / χ^β | **DONE** |
| [`fIIIa/`](fIIIa/) | \(f^{\mathrm{III}_a}\) / χ^γ | **DONE** |
| [`fIIIb/`](fIIIb/) | \(f^{\mathrm{III}_b}\) / χ^δ | **DONE** |
| [`GI/`](GI/) | \(\Gamma^{\mathrm{I}}\) / χ^ε | **DONE** |
| [`GII/`](GII/) | \(\Gamma^{\mathrm{II}}\) / χ^ζ | **DONE** (m ≡ DIRECT ≡ Path B ≡ ethS) |
| [`GIIIa/`](GIIIa/) | \(\Gamma^{\mathrm{III}_a}\) / χ^η | **DONE** (m ≡ DIRECT ≡ Path B; ladder ≡ DGEMM) |
| — | \(\Gamma^{\mathrm{III}_b}\) / χ^η RC | **DONE**\* — see [../factored_GIIIb/NOTES.md](../factored_GIIIb/NOTES.md) (pack not yet in this tree) |
| [`GIIIc/`](GIIIc/) | \(\Gamma^{\mathrm{III}_c}\) / χ^θ | **DONE** |
| [`GIVa/`](GIVa/) | \(\Gamma^{\mathrm{IV}_a}\) / χ^κ | **DONE** |
| [`GIVb/`](GIVb/) | \(\Gamma^{\mathrm{IV}_b}\) / χ^ι | **DONE**\* (fold gold any λ; Fac RC λ=0; rectangular Fac λ≠0 open) |
| [`GIVc/`](GIVc/) | \(\Gamma^{\mathrm{IV}_c}\) / χ^λ | **DONE** |

**Master table:** [../DIAGRAM_LOCK_STATUS.md](../DIAGRAM_LOCK_STATUS.md).

**Open (does not block gold):** AMC Pandya→CC ≡ code RC; rectangular Factorized RC at λ≠0 for IV_b.

## Run all gold smoke (emax=1, λ=2)

```bash
cd /Users/wolf/work/imsrg
./learn/amc_tts/tensor_pro_final/run_all_gold.sh
```

Requires `PYTHONPATH=build` and a built `pyIMSRG`.

## Packaging rule (every compare)

1. Name reduced vs unreduced on **both** sides ([REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md)).
2. m-scheme is always physical.
3. Ratios \(\sim\hat J^{\pm1}\) → packaging bug, not a wrong 6j.

## Γ^{IV_a} highlight

One Path B for **any λ**. Reduced-tensor transpose in VI_II:

\[
\bar\chi^{J_0 J_1}
=
h_\Omega\,(-1)^{J_0+J_1}
\sum\mathrm{occ}\,
\bar\Omega^{J_1 J_0}\,
\bar\Gamma^{J_1}.
\]

Details: [`GIVa/03_pathB/EQUATION.md`](GIVa/03_pathB/EQUATION.md), [../LESSONS.md](../LESSONS.md) §Reduced tensor transpose.
