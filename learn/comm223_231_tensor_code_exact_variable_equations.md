# comm223_231 Tensor: Code-Exact Variable Equations

This file is intentionally code-faithful to:
- [src/FactorizedDoubleCommutator_eths.cc](src/FactorizedDoubleCommutator_eths.cc)

It uses the exact variable names from code and the exact argument ordering used in matrix-element calls.

## Scope
- Function: comm223_231_chi1b_tensor(Eta, Gamma, Z)
- Function: comm223_231_chi2b_tensor(Eta, Gamma, Z)
- Tensor branch condition from comm223_231:
  - Eta.GetJRank() != 0
  - Gamma.GetJRank() == 0
  - Z.GetJRank() == 0

## Conventions (as coded)
- lambda = Eta.GetJRank()
- hat_lambda_inv = 1 / sqrt(2*lambda + 1)
- phase(x) means Z.modelspace->phase(x)

## A) comm223_231_chi1b_tensor

### A1) Chi_alpha(d,e)
Loop variables in code:
- outer: d, e
- inner: a, b, c, J0, J1

Code-equation:

$$
\mathrm{Chi\_alpha}(d,e)
=
\frac{1}{od.j2+1}
\sum_{a,b,c,J0,J1}
\frac{1}{2}
\,\mathrm{occfactor}(a,b,c,d)
\,\mathrm{phase}(J0+J1+lambda)
\,\mathrm{hat\_lambda\_inv}
\,\mathrm{Eta.TwoBody.GetTBME\_J}(J0,J1,c,d,a,b)
\,\mathrm{Eta.TwoBody.GetTBME\_J}(J1,J0,a,b,c,e)
$$

with

$$
\mathrm{occfactor}(a,b,c,d)
=
(1-n_a)(1-n_b)n_c n_d - n_a n_b (1-n_c)(1-n_d)
$$

and triangle filter:

$$
\mathrm{AngMom::Triangle}(J0,J1,lambda)
$$

Symmetry fill in code:

$$
\mathrm{if}\ d\neq e:\quad \mathrm{Chi\_alpha}(e,d) \mathrel{+}= \mathrm{Chi\_alpha}(d,e)
$$

### A2) Chi_beta(d,e)
Loop variables in code:
- outer: d, e
- inner: a, b, c, J0, J1

Code-equation:

$$
\mathrm{Chi\_beta}(d,e)
=
\sum_{a,b,c,J0,J1}
\frac{1}{2}
\,\mathrm{pref\_phase}(e)
\,\mathrm{occfactor}(a,b,c,d)
\,\mathrm{phase}(J0 + jc/2)
\,\sqrt{(2J0+1)(2J1+1)}
\,\mathrm{SixJ}(J0,J1,lambda,je/2,jd/2,jc/2)
$$

$$
\times\mathrm{Eta.TwoBody.GetTBME\_J}(J0,J1,c,d,a,b)
\times\mathrm{Gamma.TwoBody.GetTBME\_J}(J1,J1,a,b,c,e)
$$

where

$$
\mathrm{pref\_phase}(e) = \mathrm{phase}(je/2 + lambda)
$$

### A3) zI and zII at output indices (p,q)

Loop variables in code:
- outer: p, q
- inner: a, b
- zI sum: J0
- zII sum: J0_II, J1

#### zI

$$
zI
=
\sum_{a,b,J0}
(2J0+1)
\,\mathrm{Chi\_alpha}(a,b)
\left[
\mathrm{Gamma.TwoBody.GetTBME\_J}(J0,J0,b,p,a,q)
+
\mathrm{Gamma.TwoBody.GetTBME\_J}(J0,J0,a,p,b,q)
\right]
$$

#### zII

$$
zII
=
\sum_{a,b,J0\_II,J1}
\sqrt{(2J0\_II+1)(2J1+1)}\,\mathrm{hat\_lambda\_inv}
\,\mathrm{phase}(J1 + ja/2)
\,\mathrm{SixJ}(J1,lambda,J0\_II,jb/2,jp/2,ja/2)
$$

$$
\times\mathrm{Chi\_beta}(a,b)
\times\mathrm{Eta.TwoBody.GetTBME\_J}(J0\_II,J1,b,p,a,q)
$$

$$
-
\sum_{a,b,J0\_II,J1}
\sqrt{(2J0\_II+1)(2J1+1)}\,\mathrm{hat\_lambda\_inv}
\,\mathrm{phase}(J1 + jb/2)
\,\mathrm{SixJ}(J1,lambda,J0\_II,ja/2,jp/2,jb/2)
$$

$$
\times\mathrm{Chi\_beta}(b,a)
\times\mathrm{Eta.TwoBody.GetTBME\_J}(J0\_II,J1,a,p,b,q)
$$

#### Final update in code

$$
\mathrm{pref\_I}=\frac{1}{op.j2+1},\qquad
\mathrm{pref\_{II}}=\frac{\mathrm{phase}(op.j2/2)}{op.j2+1}
$$

$$
Z.OneBody(p,q) \mathrel{+}= \mathrm{pref\_I}\,zI + \mathrm{pref\_{II}}\,zII
$$

$$
\mathrm{if}\ p\neq q:\quad
Z.OneBody(q,p) \mathrel{+}= hZ\,\left(\mathrm{pref\_I}\,zI + \mathrm{pref\_{II}}\,zII\right)
$$

## B) comm223_231_chi2b_tensor

### B1) barred_tbme(Op,a,j,k,b,J0,J1)

Code uses:

$$
\mathrm{barred\_tbme}
=
\mathrm{pref\_phase}\cdot\mathrm{pref\_hat}
\sum_{J2,J3,j0}
\mathrm{phase}(J2)
\,\sqrt{(2J2+1)(2J3+1)}(2j0+1)
\,\mathrm{SixJ}(lam,J3,J2,ja/2,jj/2,j0)
$$

$$
\times\mathrm{SixJ}(ja/2,jb/2,J0,jk/2,j0,J3)
\times\mathrm{SixJ}(J1,J0,lam,j0,jj/2,jk/2)
	imes\mathrm{Op.TwoBody.GetTBME\_J}(J2,J3,a,b,j,k)
$$

where

$$
\mathrm{pref\_phase}=\mathrm{phase}(J0 + (ja+jk)/2 + lam),
\quad
\mathrm{pref\_hat}=\sqrt{(2J0+1)(2J1+1)}
$$

Code note: barred indices are interpreted in cross-coupled order, so
`(a,j,k,b)` in the helper maps to underlying TBME legs `(a,b;j,k)`.

### B2) bar_chi_gamma(i,l,k,j,J0)

$$
\mathrm{bar\_chi\_gamma}
=
\frac{\mathrm{phase}(J0)}{2J0+1}
\sum_{a,b,J2}
\mathrm{occ}\_{\gamma}(a,b,l,k)
\,\mathrm{phase}(J2+lambda)
\,\mathrm{hat\_lambda\_inv}
\,\mathrm{barred\_tbme}(Eta,i,b,a,j,J0,J2)
\,\mathrm{barred\_tbme}(Eta,a,l,k,b,J2,J0)
$$

with

$$
\mathrm{occ}\_{\gamma}(a,b,l,k)=n_a(1-n_b)n_l(1-n_k)-(1-n_a)n_b(1-n_l)n_k
$$

### B3) chi_delta(i,j,k,l,J0)

$$
\mathrm{chi\_delta}
=
\frac{1}{4}\frac{\mathrm{phase}(J0)}{2J0+1}
\sum_{m,n,J2}
\mathrm{occ}\_{\delta}(i,j,m,n)
\,\mathrm{phase}(J2+lambda)
\,\mathrm{hat\_lambda\_inv}
\,\mathrm{Eta.TwoBody.GetTBME\_J}(J0,J2,i,j,m,n)
\,\mathrm{Eta.TwoBody.GetTBME\_J}(J2,J0,m,n,k,l)
$$

with

$$
\mathrm{occ}\_{\delta}(i,j,m,n)=(1-n_i)(1-n_j)n_m n_n - n_i n_j(1-n_m)(1-n_n)
$$

### B4) Final zIIIa and zIIIb at (p,q)

$$
zIIIa
=
\sum_{a,b,c,J0}(2J0+1)
\left[
\mathrm{bar\_chi\_gamma}(p,c,a,b,J0)\cdot\mathrm{barred\_tbme}(Gamma,a,b,q,c,J0,J0)
-
\mathrm{bar\_chi\_gamma}(c,q,a,b,J0)\cdot\mathrm{barred\_tbme}(Gamma,a,b,c,p,J0,J0)
\right]
$$

$$
zIIIb
=
\sum_{a,b,c,J0}(2J0+1)
\left[
\mathrm{chi\_delta}(c,p,a,b,J0)\cdot\mathrm{Gamma.TwoBody.GetTBME\_J}(J0,J0,a,b,c,q)
-
\mathrm{Gamma.TwoBody.GetTBME\_J}(J0,J0,c,p,a,b)\cdot\mathrm{chi\_delta}(a,b,c,q,J0)
\right]
$$

Final one-body update:

$$
\mathrm{pref}=\frac{1}{op.j2+1}
$$

$$
Z.OneBody(p,q) \mathrel{+}= \mathrm{pref}\,(zIIIa+zIIIb)
$$

$$
\mathrm{if}\ p\neq q:\quad Z.OneBody(q,p) \mathrel{+}= hZ\,\mathrm{pref}\,(zIIIa+zIIIb)
$$

## C) Direct Code Anchors
- chi1b tensor function starts at [src/FactorizedDoubleCommutator_eths.cc](src/FactorizedDoubleCommutator_eths.cc#L114)
- chi2b tensor function starts at [src/FactorizedDoubleCommutator_eths.cc](src/FactorizedDoubleCommutator_eths.cc#L360)
