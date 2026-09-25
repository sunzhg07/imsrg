# He8 paper EOM (emax=3, EM1.8/2.0, ħω=16, NO2B)

Hs names are built from the calculation parameters and stored in `hs_cache/`:

```
Hs_{nucleus}_{LECs}_hw{hw}_emax{emax}.dat              # closed-shell SR (val = nucleus)
Hs_{nucleus}_{val}_{LECs}_hw{hw}_emax{emax}.dat        # VS / MR (val ≠ nucleus)
```

Examples:

```
hs_cache/Hs_He8_EM1.8_2.0_hw16_emax3.dat
hs_cache/Hs_He8_pshell_EM1.8_2.0_hw16_emax3.dat
```

EOM JSON goes to `eom_result/{nucleus}_{LECs}_hw{hw}_emax{emax}_{job}.json`.

RDMs are still `run2/he8.ref` and `run2/he8/*.ref`. Override LECs with `LECS=...`.

```bash
cd run_paper
export PYTHONPATH=../build

python3 he8_sr_eom.py                 # 0+  (PHASE=1 rebuilds Hs)
JRANK=2 python3 he8_sr_eom.py         # 2+
python3 he8_sr_eom.py JOB=fci         # SR Hs → UndoNO → FCI .snt (no EOM)

python3 he8_mr_eom.py JOB=mr0gs
python3 he8_mr_eom.py JOB=mr0x
python3 he8_mr_eom.py JOB=mr2gs
python3 he8_mr_eom.py JOB=mr2x
MAX_ITER=100 python3 he8_mr_eom.py JOB=t02
python3 he8_mr_eom.py JOB=t20

# VS Hs → packaging (no EOM); both writers share the same post-decouple Hs
python3 he8_hs_fci.py                 # HF + VS (He8, p-shell) → UndoNO → KSHELL FCI .snt
python3 he8_mr_eom.py JOB=pshell      # He4-core p-shell .snt (KSHELL)
python3 he8_mr_eom.py JOB=fci         # vacuum FCI .snt (same as he8_hs_fci.py)
python3 he8_mr_eom.py JOB=snt         # both from one decouple

python3 he8_spectrum.py

# leftover / EOM benches (emax=2 He8 / He4+p-shell). Grid via LAMBDAS, TS, PARITIES.
EMAX=2 LAMBDAS=0,1,2,3 TS=0,1,2 PARITIES=0,1 python3 test_computenorm_tensor.py
EMAX=2 LAMBDAS=0,1,2,3 TS=0,1,2 PARITIES=0,1 python3 test_proj_preserves_norm_tensor.py
EMAX=2 LAMBDAS=0,1,2,3 TS=0,1,2 PARITIES=0,1 python3 test_hab_hba_tensor.py
python3 test_arnoldi_tensor.py              # ArnoldiSolveH2 Rayleigh = ⟨H⟩
python3 test_fact_unfact_223_tts.py         # nested 223st→132+231+232 vs ethS
LAMBDAS=3 TS=1 PARITIES=0,1 python3 test_computenorm_tensor.py
```
