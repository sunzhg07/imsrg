# AMC Yutsis / nj reduction — how to show χ^γ path equivalence

## How AMC reduces CG → 6j/9j/12j

AMC does **not** take 6j/9j as input. Flow:

1. m-scheme equation → Clebsch–Gordan network  
2. **Yutsis graph** reduction (`amc/yutsis/`) → 6j (and sometimes 3j constraints)  
3. Optional post-process: `--collect-ninejs` merges **3×6j → 1×9j** when an auxiliary \(j\) has phase \(2x\) and factor \((2x+1)\)  
4. Internal code also has `collect_twelvejfirsts` (12j of the first kind) — **not** exposed on the CLI  

Docs: `amc/docs/ug.rst` (`--collect-ninejs`, `--keep-threejs`, `--wet-convention`).

```bash
cd /Users/wolf/work/amc
PYTHONPATH=. python3 -m amc -o out.tex input.txt
PYTHONPATH=. python3 -m amc --collect-ninejs -o out_9j.tex input.txt
PYTHONPATH=. python3 -m amc -v -o out.tex input.txt   # prints CG / Yutsis debug
```

## Equivalence inputs (this folder)

| File | Meaning |
|---|---|
| `input/chi_gamma_direct.txt` | Path A: m-scheme → normal \(\chi^\gamma\) (slow path) |
| `input/chi_gamma_via_pandya.txt` | Path B: \(\Omega\to\bar\Omega\to\bar\chi\to\chi\) |

Outputs: `output/chi_gamma_{direct,via_pandya}_{plain,ninej}.tex`

## What AMC shows (same physics, different factorization)

**Path A (direct / slow):**

- plain: **5×6j** × \(\Omega\Omega\)  
- `--collect-ninejs`: **2×6j + 1×9j** × \(\Omega\Omega\)

**Path B (Pandya pipeline), three equations:**

1. \(\bar\Omega\): tensor Pandya = **3×6j** (plain) or **1×9j** (`--collect-ninejs`)  
2. \(\bar\chi^\gamma\): **RME product** in Pandya — \(\hat\lambda^{-1}\,\bar\Omega^{J_2 J_0\lambda}\bar\Omega^{J_0 J_2\lambda}\) (no extra 6j)  
3. \(\chi^\gamma\): scalar inv Pandya = **1×6j** × \(\bar\chi\)

So tensor “Pandya × Pandya × inv” is **9j + RME + 6j** (or 3×6j + RME + 6j), not “only two 6j”.  
The slow-path \(W_1\)’s two 6j are only **one leg** of Path A’s five.

AMC does **not** auto-substitute Path B into a single Ω-only formula; equivalence is by construction of the pipeline (Pandya def + product + inv) vs direct reduction of the same m-scheme. To *prove* identity algebraically one substitutes (1)+(2) into (3) and uses recoupling; numerically, Path A χ is the oracle for the composed Path B.
