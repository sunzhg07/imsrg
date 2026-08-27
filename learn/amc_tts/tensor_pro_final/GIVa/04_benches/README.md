# $\Gamma^{\mathrm{IV}_a}$ / $\chi^\kappa$ — three-way benches (m ≡ AMC direct ≡ Path B)

Status: **DONE**

Gold: m ≡ AMC ≡ Path B ≡ ethS (λ≠0). Lessons: LESSONS.md §Reduced tensor transpose.

## Scripts (from repo `run/`)

```bash
cd /Users/wolf/work/imsrg
export PYTHONPATH=build
python3 -B run/test_chi_kappa_m_vs_amc.py 1 2   # emax λ — see script docstring
python3 -B run/test_chi_kappa_pathB_vs_direct.py 1 2   # emax λ — see script docstring
python3 -B run/test_G4a_Wbra_mscheme.py 1 2   # emax λ — see script docstring
python3 -B run/test_G4a_pathB_mscheme.py 1 2   # emax λ — see script docstring
python3 -B run/test_tts_GIVa_eths_vs_pathB.py 1 2   # emax λ — see script docstring
```

Linked copies under `scripts/` point at the live `run/` files.

