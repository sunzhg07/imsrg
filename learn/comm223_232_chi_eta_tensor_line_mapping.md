# comm223_232 chi_eta tensor line mapping

Scope: only chi^eta path in comm223_232_chi2b.

## Equation symbol to code mapping

1. chi^eta branch selector
- Equation role: switch scalar vs tensor eta route.
- Code symbol: tensor_eta_case, lambda, hat_lambda_inv.
- Lines: src/FactorizedDoubleCommutator_eths.cc:1094-1097.

2. Barred tensor Omega for eta (with adcb convention)
- Equation role: barred Omega recoupling for tensor lambda.
- Code symbol: barred_eta_tensor(a,j,k,b,J0,J1).
- Direct TBME call carrying adcb interpretation: Omega_{abjk} via Eta.TwoBody.GetTBME_J(J2,J3,a,b,j,k).
- Lines: src/FactorizedDoubleCommutator_eths.cc:1101-1167.

3. chi^eta in cross-coupled space (bar chi)
- Equation role: bar(chi^eta) ~ sum occ * bar(Omega) * bar(Omega).
- Scalar fallback: barCHI_III[ch_cc] = bar_Eta[ch_cc] * nnnbar_Eta.
- Tensor route: explicit loops over ibra_cc, iket_cc, a, b, J2 and bo1/bo2 from barred_eta_tensor.
- Output: barCHI_III[ch_cc](ibra_cc, iket_cc) = phase(J_cc) * sum_eta / (2J_cc+1).
- Lines: src/FactorizedDoubleCommutator_eths.cc:1329-1396.

4. Reverse Pandya of chi^eta to Chi_III_Op, direct ijkl leg
- Equation role: inverse transform from barCHI_III to direct-channel Chi_III_Op.
- Barred access: me1 = barCHI_III[ch_cc](indx_il, indx_kj).
- Scalar route: commij -= (2Jprime+1) * sixj * me1.
- Tensor route: commij += pref * rec * me1, with rec built from three 6j symbols and j0 loop.
- Lines: src/FactorizedDoubleCommutator_eths.cc:1487-1513.

5. Reverse Pandya of chi^eta to Chi_III_Op, exchanged i<->j leg
- Equation role: exchange contribution for antisymmetrized combination.
- Barred access: me1 = barCHI_III[ch_cc](indx_lj, indx_ik).
- Scalar route: commji -= (2Jprime+1) * sixj * me1.
- Tensor route: commji += pref * rec * me1, with exchanged three-6j structure.
- Lines: src/FactorizedDoubleCommutator_eths.cc:1540-1565.

6. Final scalar object assembly (unchanged downstream)
- Equation role: build direct 2b chi^eta object then contract in scalar matrix products.
- Code: zijkl from commij/commji, then add to Chi_III_Op.GetMatrix(ch,ch).
- Lines: src/FactorizedDoubleCommutator_eths.cc:1569-1579.

## Loop-index map for tensor chi^eta core

1. Cross-coupled bra/ket loops
- ibra_cc -> (i,l)
- iket_cc -> (k,j)
- Lines: src/FactorizedDoubleCommutator_eths.cc:1333-1353.

2. Internal sums
- a,b over all_orbits, J2 over coupled range constrained by Triangle(J_cc,J2,lambda).
- Lines: src/FactorizedDoubleCommutator_eths.cc:1362-1383.

3. Two barred Omegas in product
- bo1 = barred_eta_tensor(i,b,a,j,J_cc,J2)
- bo2 = barred_eta_tensor(a,l,k,b,J2,J_cc)
- Lines: src/FactorizedDoubleCommutator_eths.cc:1384-1385.

## adcb convention checkpoints

1. Tensor barred helper comment and TBME leg choice
- Comment states barred (a,j,k,b) maps to direct legs (a,b;j,k).
- Implementation anchor: Eta.TwoBody.GetTBME_J(J2,J3,a,b,j,k).
- Lines: src/FactorizedDoubleCommutator_eths.cc:1100-1102 and 1145.

2. Reverse-Pandya barred read ordering
- First leg uses (il,kj): barCHI_III(ch_cc)(indx_il,indx_kj).
- Exchange leg uses (lj,ik): barCHI_III(ch_cc)(indx_lj,indx_ik).
- Lines: src/FactorizedDoubleCommutator_eths.cc:1487 and 1540.

## Reference equation anchor used

- diag2 compact chi^eta (barred) form:
  amc/examples/sample_output/diag2_compact.tex:80.
