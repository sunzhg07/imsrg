# $\Gamma^{\mathrm{I}}$ / $\chi^\varepsilon$ — three-way benches (m ≡ AMC direct ≡ Path B)

Status: **DONE**

Gold: m ≡ AMC ≡ Path B / DGEMM.

## Scripts (from repo `run/`)

```bash
cd /Users/wolf/work/imsrg
export PYTHONPATH=build
python3 -B run/test_chi_epsilon_mscheme.py 1 2   # emax λ — see script docstring
python3 -B run/test_tts_GI_mscheme.py 1 2   # emax λ — see script docstring
python3 -B run/test_tts_GI.py 1 2   # emax λ — see script docstring
python3 -B run/test_tts_GI_dgemm.py 1 2   # emax λ — see script docstring
```

Linked copies under `scripts/` point at the live `run/` files.

