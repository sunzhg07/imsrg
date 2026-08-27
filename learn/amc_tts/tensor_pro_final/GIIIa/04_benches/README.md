# $\Gamma^{\mathrm{III}_a}$ / $\chi^\eta$ — three-way benches (m ≡ AMC direct ≡ Path B)

Status: **DONE** — χ m ≡ DIRECT ≡ Path B; ladder m ≡ from_chi ≡ DGEMM

(Optional TTS dual-oracle only; not required for gold lock.)

## Scripts (from repo `run/`)

```bash
cd /Users/wolf/work/imsrg
export PYTHONPATH=build
python3 -B run/test_chi_eta_mscheme.py 1 2   # emax λ — see script docstring
python3 -B run/test_GIIIa_ladder_mscheme.py 1 2   # emax λ — see script docstring
python3 -B run/test_tts_GIIIa.py 1 2   # emax λ — see script docstring
```

Linked copies under `scripts/` point at the live `run/` files.

