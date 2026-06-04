# Tensor-coupled AMC setup for the factorized double commutators

## Goal

This note rewrites the same M-scheme double-commutator equations from [learn/arxiv_eq.tex](learn/arxiv_eq.tex), but now with the following operator content:

- $\Omega$ is a tensor operator,
- $H$ or $\Gamma$ is still a scalar operator,
- the final stored results in $f$ and $\Gamma$ are always scalar operators.

> **Important clarification for the actual IMSRG use case in this repository:** we do **not** want a tensor final
> one-body operator $f$ here.  The physical/output equations of interest remain **scalar**.  In particular, for the
> first factorized term one should keep the same scalar closure logic as in the arXiv equations:
> $$
> \chi^\alpha \sim [\Omega\times\Omega]^0 \quad\text{(scalar)},
> \qquad
> f^I \sim \chi^\alpha\,\Gamma \quad\text{(scalar $\times$ scalar $\to$ scalar)}.
> $$
> More generally, in the scalar-output organization one has
> - $\Omega\Omega \to$ scalar intermediate,
> - $\Gamma\Omega \to$ tensor intermediate,
> - tensor $\times$ tensor must then be coupled to rank $0$ to recover a scalar final $f$.
>
> So the correct classification to use for the production equations is:
> intermediates may be scalar or tensor, but the stored/output $f$ and $\Gamma$ objects are scalar.

Because $H$ is scalar, coupling a tensor operator $\Omega^{(\lambda)}$ with $H^{(0)}$ can generate tensor
intermediates in the reduction, but in the present application the final stored operators are projected back to
rank $0$,

$$
[\Omega^{(\lambda)}, H^{(0)}] \;\leadsto\; \text{tensor intermediates and/or scalar intermediates, then scalar final output.}
$$

The uncoupled M-scheme algebra is unchanged. What changes is the angular-momentum reduction: in J-scheme, the reduced matrix elements now have different angular momentum labels on bra and ket,

$$
\Omega_{abcd}^{J_{ab} J_{cd} \lambda},
$$

$$
C_{ijkl}^{J_{ij} J_{kl} \lambda},
$$

instead of the scalar pattern

$$
A_{abcd}^{J J 0}.
$$

AMC can handle this organization when `scalar=false` is used for $\Omega$, while the final declared output
operators $f$ and $\Gamma$ remain scalar.

## AMC declarations

A minimal AMC declaration block is

```none
declare Omega {
    mode=4,
    scalar=false,
    latex="\Omega",
}

declare Gamma {
    mode=4,
    scalar=true,
    latex="\Gamma",
}

declare n {
    mode=2,
    diagonal=true,
    latex="n",
}

declare nbar {
    mode=2,
    diagonal=true,
    latex="\bar{n}",
}

declare f1 { mode=2, scalar=true, latex="f^{(I)}" }
declare f2 { mode=2, scalar=true, latex="f^{(II)}" }
declare f3a { mode=2, scalar=true, latex="f^{(III_a)}" }
declare f3b { mode=2, scalar=true, latex="f^{(III_b)}" }

declare G1 { mode=4, scalar=true, latex="\Gamma^{(I)}" }
declare G2 { mode=4, scalar=true, latex="\Gamma^{(II)}" }
declare G3a { mode=4, scalar=true, latex="\Gamma^{(III_a)}" }
declare G3b { mode=4, scalar=true, latex="\Gamma^{(III_b)}" }
declare G3c { mode=4, scalar=true, latex="\Gamma^{(III_c)}" }
declare G4a { mode=4, scalar=true, latex="\Gamma^{(IV_a)}" }
declare G4b { mode=4, scalar=true, latex="\Gamma^{(IV_b)}" }
declare G4c { mode=4, scalar=true, latex="\Gamma^{(IV_c)}" }
```

Here `Gamma` is scalar, `Omega` is tensor, and the stored outputs are scalar operators.  Tensor character may appear
in intermediate coupled objects, but the final stored one-body output is scalar, schematically of the form

$$
f_{ij}^{0},
$$

and the stored two-body output is likewise scalar,

$$
\Gamma_{ijkl}^{J J 0}.
$$

## One-body M-scheme equations

These equations are identical in uncoupled form to the scalar case.

### $f^{I}_{ij}$

$$
f^{I}_{ij} = \frac{1}{2} \sum_{abcde}
\left( \bar n_a \bar n_b n_c n_d - n_a n_b \bar n_c \bar n_d \right)
\left(
\Omega_{cdab} \Omega_{abce} \Gamma_{eidj}
+
\Omega_{cdab} \Omega_{abce} \Gamma_{diej}
\right).
$$

AMC input:

```none
f1_ij = 1/2 * sum_abcde(
    (nbar_a*nbar_b*n_c*n_d - n_a*n_b*nbar_c*nbar_d)
    *
    (Omega_cdab*Omega_abce*Gamma_eidj + Omega_cdab*Omega_abce*Gamma_diej)
);
```

### $f^{II}_{ij}$

$$
f^{II}_{ij} = \frac{1}{2} \sum_{abcde}
\left( n_a n_b \bar n_c \bar n_e - \bar n_a \bar n_b n_c n_e \right)
\left(
\Gamma_{cdab} \Omega_{abce} \Omega_{eidj}
-
\Gamma_{cdab} \Omega_{abce} \Omega_{diej}
\right).
$$

AMC input:

```none
f2_ij = 1/2 * sum_abcde(
    (n_a*n_b*nbar_c*nbar_e - nbar_a*nbar_b*n_c*n_e)
    *
    (Gamma_cdab*Omega_abce*Omega_eidj - Gamma_cdab*Omega_abce*Omega_diej)
);
```

### $f^{III_a}_{ij}$

$$
f^{III_a}_{ij} = \sum_{abcde}
\left( \bar n_a \bar n_b n_c n_d - n_a n_b \bar n_c \bar n_d \right)
\left(
\Omega_{abcd} \Omega_{idae} \Gamma_{cejb}
-
\Omega_{abcd} \Omega_{edaj} \Gamma_{cieb}
\right).
$$

AMC input:

```none
f3a_ij = sum_abcde(
    (nbar_a*nbar_b*n_c*n_d - n_a*n_b*nbar_c*nbar_d)
    *
    (Omega_abcd*Omega_idae*Gamma_cejb - Omega_abcd*Omega_edaj*Gamma_cieb)
);
```

### $f^{III_b}_{ij}$

$$
f^{III_b}_{ij} = \frac{1}{4} \sum_{abcde}
\left( \bar n_a \bar n_b n_c n_d - n_a n_b \bar n_c \bar n_d \right)
\left(
\Omega_{abcd} \Omega_{cdej} \Gamma_{eiab}
-
\Omega_{abcd} \Omega_{eiab} \Gamma_{cdej}
\right).
$$

AMC input:

```none
f3b_ij = 1/4 * sum_abcde(
    (nbar_a*nbar_b*n_c*n_d - n_a*n_b*nbar_c*nbar_d)
    *
    (Omega_abcd*Omega_cdej*Gamma_eiab - Omega_abcd*Omega_eiab*Gamma_cdej)
);
```

## Two-body M-scheme equations

Again, the uncoupled equations are unchanged. Only the reduced J-coupled structure changes because $\Omega$ is now a tensor operator.

### $\Gamma^{I}_{ijkl}$

$$
\Gamma^{I}_{ijkl} = \frac{1}{2} \sum_{abcd}
\left( \bar n_a \bar n_b n_c + n_a n_b \bar n_c \right)
\left\{
(1-\hat P_{ij}) \Omega_{ciab} \Omega_{abcd} \Gamma_{djkl}
+
(1-\hat P_{kl}) \Omega_{cdab} \Omega_{abcl} \Gamma_{ijkd}
\right\}.
$$

AMC input:

```none
G1_ijkl = 1/2 * sum_abcd(
    (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
    *
    (P(i/j)*Omega_ciab*Omega_abcd*Gamma_djkl + P(k/l)*Omega_cdab*Omega_abcl*Gamma_ijkd)
);
```

### $\Gamma^{II}_{ijkl}$

$$
\Gamma^{II}_{ijkl} = -\frac{1}{2} \sum_{abcd}
\left( \bar n_a \bar n_b n_c + \bar n_c n_a n_b \right)
\left\{
(1-\hat P_{ij}) \Omega_{cjab} \Gamma_{abcd} \Omega_{idkl}
+
(1-\hat P_{kl}) \Gamma_{cdab} \Omega_{abcl} \Omega_{ijkd}
\right\}.
$$

AMC input:

```none
G2_ijkl = -1/2 * sum_abcd(
    (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
    *
    (P(i/j)*Omega_cjab*Gamma_abcd*Omega_idkl + P(k/l)*Gamma_cdab*Omega_abcl*Omega_ijkd)
);
```

### $\Gamma^{III_a}_{ijkl}$

$$
\Gamma^{III_a}_{ijkl} = - \sum_{abcd}
\left( \bar n_c \bar n_d n_a + n_a \bar n_c \bar n_d \right)
\left\{
(1-\hat P_{ij}) \Omega_{ajcd} \Omega_{idab} \Gamma_{cbkl}
+
(1-\hat P_{kl}) \Omega_{cdka} \Omega_{bacl} \Gamma_{ijbd}
\right\}.
$$

AMC input:

```none
G3a_ijkl = -1 * sum_abcd(
    (nbar_c*nbar_d*n_a + n_a*nbar_c*nbar_d)
    *
    (P(i/j)*Omega_ajcd*Omega_idab*Gamma_cbkl + P(k/l)*Omega_cdka*Omega_bacl*Gamma_ijbd)
);
```

### $\Gamma^{III_b}_{ijkl}$

$$
\Gamma^{III_b}_{ijkl} = - \sum_{abcd}
\left( \bar n_b n_c n_d + n_b \bar n_c \bar n_d \right)
(1-\hat P_{ij})(1-\hat P_{kl})
\left(
\Omega_{dcbk} \Omega_{biac} \Gamma_{jald}
+
\Omega_{jcbd} \Omega_{balc} \Gamma_{diak}
\right).
$$

AMC input:

```none
G3b_ijkl = -1 * sum_abcd(
    (nbar_b*n_c*n_d + n_b*nbar_c*nbar_d)
    * P(i/j) * P(k/l)
    *
    (Omega_dcbk*Omega_biac*Gamma_jald + Omega_jcbd*Omega_balc*Gamma_diak)
);
```

### $\Gamma^{III_c}_{ijkl}$

$$
\Gamma^{III_c}_{ijkl} = -\frac{1}{2} \sum_{abcd}
\left( \bar n_a \bar n_b n_c + n_a n_b \bar n_c \right)
(1-\hat P_{ij})(1-\hat P_{kl})
\left(
\Omega_{abcl} \Omega_{idab} \Gamma_{cjkd}
+
\Omega_{icab} \Omega_{abdl} \Gamma_{djkc}
\right).
$$

AMC input:

```none
G3c_ijkl = -1/2 * sum_abcd(
    (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
    * P(i/j) * P(k/l)
    *
    (Omega_abcl*Omega_idab*Gamma_cjkd + Omega_icab*Omega_abdl*Gamma_djkc)
);
```

### $\Gamma^{IV_a}_{ijkl}$

$$
\Gamma^{IV_a}_{ijkl} = - \sum_{abcd}
\left( \bar n_c \bar n_d n_a + n_c n_d \bar n_a \right)
\left\{
(1-\hat P_{ij}) \Omega_{aicd} \Omega_{dbkl} \Gamma_{jcba}
+
(1-\hat P_{kl}) \Omega_{dcak} \Omega_{ijcb} \Gamma_{bald}
\right\}.
$$

AMC input:

```none
G4a_ijkl = -1 * sum_abcd(
    (nbar_c*nbar_d*n_a + n_c*n_d*nbar_a)
    *
    (P(i/j)*Omega_aicd*Omega_dbkl*Gamma_jcba + P(k/l)*Omega_dcak*Omega_ijcb*Gamma_bald)
);
```

### $\Gamma^{IV_b}_{ijkl}$

$$
\Gamma^{IV_b}_{ijkl} = (1-\hat P_{ij})(1-\hat P_{kl})
\sum_{abcd}
\left( \bar n_a n_b \bar n_c + n_a \bar n_b n_c \right)
\left(
\Omega_{bica} \Omega_{jcld} \Gamma_{dabk}
+
\Omega_{cabk} \Omega_{jdlc} \Gamma_{bida}
\right).
$$

AMC input:

```none
G4b_ijkl = sum_abcd(
    (nbar_a*n_b*nbar_c + n_a*nbar_b*n_c)
    * P(i/j) * P(k/l)
    *
    (Omega_bica*Omega_jcld*Gamma_dabk + Omega_cabk*Omega_jdlc*Gamma_bida)
);
```

### $\Gamma^{IV_c}_{ijkl}$

$$
\Gamma^{IV_c}_{ijkl} = \frac{1}{2} (1-\hat P_{ij})(1-\hat P_{kl})
\sum_{abcd}
\left( \bar n_a \bar n_b n_d + n_a n_b \bar n_d \right)
\left(
\Omega_{abld} \Omega_{djck} \Gamma_{icab}
+
\Omega_{idab} \Omega_{cjdk} \Gamma_{ablc}
\right).
$$

AMC input:

```none
G4c_ijkl = 1/2 * sum_abcd(
    (nbar_a*nbar_b*n_d + n_a*n_b*nbar_d)
    * P(i/j) * P(k/l)
    *
    (Omega_abld*Omega_djck*Gamma_icab + Omega_idab*Omega_cjdk*Gamma_ablc)
);
```

## Factorized scheme with arXiv-style intermediates

To follow the organization in [learn/arxiv_eq.tex](learn/arxiv_eq.tex), it is useful to keep the same two classes of intermediates:

1. products of two $\Omega$ operators,
2. products of one scalar $H$ or $\Gamma$ with one tensor $\Omega$.

With $\Omega$ nonscalar and $\Gamma$ scalar, the natural tensor character of these intermediates is

$$
\Omega \times \Omega \longrightarrow \text{scalar intermediate},
$$

$$
\Gamma \times \Omega \longrightarrow \text{tensor intermediate of rank } \lambda.
$$

So the organization we want to keep is:

- scalar intermediates built from $\Omega\Omega$,
- tensor intermediates built from $\Gamma\Omega$,
- the stored final operators are scalar, even if AMC exposes tensor-coupled intermediate objects along the way.

### Scalar intermediates from $\Omega\Omega$

The direct tensor-coupled analogues of the scalar arXiv intermediates are

$$
\chi^{\alpha}_{ij}
= \sum_{abcJ}
\frac{\hat J^2}{2\hat j_i^2}
\left(
\bar n_a \bar n_b n_c n_i - n_a n_b \bar n_c \bar n_i
\right)
\left[ \Omega_{ciab}^{J J \lambda} \times \Omega_{abcj}^{J J \lambda} \right]^0,
$$

$$
\chi^{\delta J}_{ijkl}
= \sum_{ab}
\frac{\hat J^2}{4}
\left(
n_a n_b \bar n_k \bar n_l - \bar n_a \bar n_b n_k n_l
\right)
\left[ \Omega_{ijab}^{J J \lambda} \times \Omega_{abkl}^{J J \lambda} \right]^0,
$$

$$
\chi^{\epsilon}_{ij}
= \frac{1}{2\hat j_j^2}
\sum_{abcJ}
\hat J^2
\left(
\bar n_a \bar n_b n_c + n_a n_b \bar n_c
\right)
\left[ \Omega_{ciab}^{J J \lambda} \times \Omega_{abcj}^{J J \lambda} \right]^0,
$$

$$
\bar\chi^{\eta J}_{i\bar j k\bar l}
= \sum_{ab}
\left(
\bar n_a n_b \bar n_k + n_a \bar n_b n_k
\right)
\left[ \bar\Omega^{J}_{i\bar j a\bar b} \times \bar\Omega^{J}_{a\bar b k\bar l} \right]^0,
$$

$$
\bar\chi^{\theta J}_{i\bar l k\bar j}
= \sum_{ab}
\left(
n_a n_b \bar n_k + \bar n_a \bar n_b n_k + n_a n_b \bar n_j + \bar n_a \bar n_b n_j
\right)
\left[ \bar\Omega^{J}_{i\bar l a\bar b} \times \bar\Omega^{J}_{a\bar b k\bar j} \right]^0.
$$

These are scalar with respect to the external tensor rank. They carry only the recoupling labels needed to connect the two $\Omega$ factors.

### Tensor intermediates from $\Gamma\Omega$

The tensor analogues of the mixed intermediates are

$$
\Xi^{\beta,\lambda}_{ij}
= \sum_{abcJ}
\frac{\hat J^2}{2\hat j_i^2}
\left(
\bar n_a \bar n_b n_c n_i - n_a n_b \bar n_c \bar n_i
\right)
\Gamma_{ciab}^{J}
\Omega_{abcj}^{J J \lambda},
$$

$$
\Xi^{\zeta,\lambda}_{ij}
= \frac{1}{2\hat j_j^2}
\sum_{abcJ}
\hat J^2
\left(
\bar n_a \bar n_b n_c + n_a n_b \bar n_c
\right)
\Gamma_{aibc}^{J}
\Omega_{bcaj}^{J J \lambda},
$$

$$
\bar\Xi^{\iota,\lambda J}_{i\bar j k\bar l}
= \sum_{ab}
\left(
\bar n_a n_b \bar n_k + n_a \bar n_b n_k
\right)
\bar\Gamma^{J}_{i\bar j a\bar b}
\bar\Omega^{J}_{a\bar b k\bar l},
$$

$$
\overline{\overline{\Xi}}^{\kappa,\lambda J}_{i\bar j k\bar l}
= \sum_{ab}
\left(
\bar n_a n_b n_l + n_a \bar n_b \bar n_l
\right)
\overline{\overline\Gamma}^{J}_{b\bar a k\bar l}
\overline{\overline\Omega}^{J}_{i\bar j b\bar a},
$$

$$
\bar\Xi^{\lambda,\lambda J}_{i\bar l k\bar j}
= \sum_{ab}
\left(
\bar n_a \bar n_b n_l + n_a n_b \bar n_l
\right)
\Gamma^{J}_{ijab}
\bar\Omega^{J}_{a\bar b k\bar j}
+ \sum_{ab}
\left(
\bar n_a \bar n_b n_j + n_a n_b \bar n_j
\right)
\bar\Omega^{J}_{i\bar l a\bar b}
\Gamma^{J}_{abkl}.
$$

These objects carry the tensor rank $\lambda$ of $\Omega$ at the intermediate stage only.

### One-body factorization pattern and interpretation

With these definitions, the one-body equations can be organized in the same spirit as the arXiv equations.

**Important:** the displayed LaTeX equations below are being kept as AMC-style recoupling templates / copied-output-style
equations.  They are **not** to be manually rewritten here.  Our scalar-output requirement should therefore be read as
an interpretation of how these tensor objects are used in the final storage/closure, not as a change to the displayed
AMC equations themselves.

For storage/physics interpretation only:

- term I: scalar $\chi^\alpha$ times scalar $\Gamma$,
- term II: tensor $\Xi^\beta$ times tensor $\Omega$, coupled to total rank $0$,
- similarly for the remaining terms.

So the $f^{I,\lambda}$, $f^{II,\lambda}$, etc. shown below should be understood as AMC recoupling templates rather than
as a license to edit the AMC LaTeX into a different displayed equation.

$$
f^{I,\lambda}_{ij}
= \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abJ}
\hat J^2
\chi^{\alpha}_{ab}
\Omega^{J J \lambda}_{biaj},
$$

$$
f^{II,\lambda}_{ij}
= \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abJ}
\hat J^2
\left(
\Xi^{\beta,\lambda}_{ab} - (-1)^{j_a+j_b-J} \Xi^{\beta,\lambda}_{ba}
\right),
$$

$$
f^{III_a,\lambda}_{ij}
= \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abcJ}
\left(
\bar\chi^{\eta J}_{i\bar c a\bar b} \; \bar\Xi^{\gamma,\lambda J}_{a\bar b j\bar c}
-
\bar\chi^{\eta J}_{c\bar j a\bar b} \; \bar\Xi^{\gamma,\lambda J}_{a\bar b c\bar i}
\right),
$$

$$
f^{III_b,\lambda}_{ij}
= \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abcJ}
\left(
\chi^{\delta J}_{ciab} \; \Omega^{J J \lambda}_{abcj}
-
\chi^{\delta J}_{abcj} \; \Omega^{J J \lambda}_{ciab}
\right).
$$

The important point is not the exact phase convention in these schematic equations, but the separation:

- the $\Omega\Omega$ block is scalar,
- the block carrying one $\Omega$ and one $\Gamma$ is tensor,
- the final stored result is scalar, even if the displayed AMC-style formula carries tensor labels at the intermediate/template level.

Equivalently: whenever a mixed tensor intermediate appears, the remaining tensor factor must be coupled/projected with
it to total rank $0$ in the final storage logic.

In particular, the first term should remain a scalar closure of the form

$$
f^I_{ij} \sim \sum_{abJ} \chi^\alpha_{ab}\,\Gamma^J_{biaj},
$$

with both $\chi^\alpha$ and $\Gamma$ scalar, whereas mixed terms such as term II use a tensor intermediate from
$\Gamma\Omega$ which must be contracted with the remaining tensor $\Omega$ to total rank $0$.

### Two-body factorization pattern and interpretation

The same pattern extends to the two-body equations.  Again, the displayed equations below are left in AMC/template
form; the scalar-output requirement applies to how the final operators are interpreted/stored, not to rewriting the
displayed AMC LaTeX by hand.

$$
\Gamma^{I,\lambda J}_{ijkl}
= \sum_a
\left[
(1-\hat P^J_{ij}) \chi^{\epsilon}_{ai} \Omega^{J J \lambda}_{ajkl}
+
(1-\hat P^J_{kl}) \chi^{\epsilon}_{ak} \Omega^{J J \lambda}_{ijal}
\right],
$$

$$
\Gamma^{II,\lambda J}_{ijkl}
= \sum_a
\left[
(1-\hat P^J_{ij}) \Xi^{\zeta,\lambda}_{aj} \Gamma^{J}_{iakl}
-
(1-\hat P^J_{kl}) \Xi^{\zeta,\lambda}_{ak} \Gamma^{J}_{ijal}
\right],
$$

$$
\Gamma^{III_a,\lambda J}_{ijkl}
= \sum_{ab}
\left[
(1-\hat P^J_{ij}) \bar\chi^{\eta J}_{ijab} \; \bar\Xi^{\gamma,\lambda J}_{abkl}
+
(1-\hat P^J_{kl}) \bar\Xi^{\gamma,\lambda J}_{ijab} \; \bar\chi^{\eta J}_{abkl}
\right],
$$

$$
\overline{\overline\Gamma}^{III_b,\lambda J}_{j\bar l k\bar i}
= (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\overline{\overline\Omega}^{J}_{j\bar l a\bar b}
\left(
\overline{\overline\chi}^{\eta J}_{a\bar b k\bar i}
+
\overline{\overline\chi}^{\eta J}_{k\bar i b\bar a}
\right),
$$

$$
\bar\Gamma^{III_c,\lambda J}_{i\bar l k\bar j}
= \frac12 (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\bar\chi^{\theta J}_{i\bar l a\bar b}
\bar\Xi^{\gamma,\lambda J}_{a\bar b k\bar j},
$$

$$
\Gamma^{IV_a,\lambda J}_{ijkl}
= -\sum_{ab}
\left[
(1-\hat P^J_{ij}) \overline{\overline\Xi}^{\kappa,\lambda J}_{ijab} \Omega^{J}_{abkl}
+
(1-\hat P^J_{kl}) \Omega^{J}_{ijab} \overline{\overline\Xi}^{\kappa,\lambda J}_{klab}
\right],
$$

$$
\overline{\overline\Gamma}^{IV_b,\lambda J}_{j\bar l k\bar i}
= (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\overline{\overline\Omega}^{J}_{j\bar l a\bar b}
\left(
\bar\Xi^{\iota,\lambda J}_{a\bar b k\bar i}
-
\bar\Xi^{\iota,\lambda J}_{k\bar i a\bar b}
\right),
$$

$$
\bar\Gamma^{IV_c,\lambda J}_{i\bar l k\bar j}
= \frac12 (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\bar\Xi^{\lambda,\lambda J}_{i\bar l a\bar b}
\bar\Omega^{J}_{a\bar b k\bar j}.
$$

These equations should be read as recoupling templates for intermediate analysis, with the LaTeX kept verbatim in the
AMC/template style. The exact phases and recoupling coefficients must still be fixed by explicit AMC reduction term by
term.

## What changes in the J-coupled reduction

When $\Omega$ is nonscalar and $\Gamma$ is scalar, AMC can produce reduced matrix elements with tensor rank $\lambda$
at the intermediate stage. Schematically,

$$
\Omega_{abcd}^{J_1 J_2 \lambda},
$$

$$
f_{ij}^{\lambda},
$$

$$
\Gamma_{ijkl}^{J_1 J_2 \lambda}.
$$

For the present application, these tensor objects are not the stored final outputs; they are recoupling templates that
must be combined/projected to recover scalar stored equations of the form

$$
A_{abcd}^{J J 0},
$$

but instead require explicit recoupling between different bra-side and ket-side angular momenta. In practice, AMC introduces additional sums over intermediate couplings and the tensor rank $\lambda$, together with $6j$ and $9j$ symbols.

This is exactly the same mechanism seen in the nonscalar Pandya example in [amc/examples/sample_output/pandya.tex](amc/examples/sample_output/pandya.tex), where the transformed tensor matrix elements take the form

$$
\tilde V_{abcd}^{J_0 J_1 \lambda}
= \sum_{J_2 J_3 j_0}
(\text{phase})
\hat J_0 \hat J_1 \hat J_2 \hat J_3 \hat j_0^2
\times 6j \times 6j \times 6j
\times V_{abcd}^{J_2 J_3 \lambda}.
$$

So for the double-commutator equations, one should expect the reduction to differ from the naive scalar formulas in two ways:

1. tensor intermediates may appear internally with rank $\lambda$,
2. recoupling coefficients connect the left and right two-body couplings because $J_{\mathrm{bra}}$ and $J_{\mathrm{ket}}$ need not be equal, before projection back to scalar output.

## Practical note

The M-scheme equations above can be used directly as AMC input once `Omega` is declared as `scalar=false`; the final
stored output operators should still be declared scalar.  But the AMC-derived LaTeX equations themselves should be kept
verbatim rather than manually rewritten.

The most important conceptual point is:

$$
\text{same uncoupled algebra} \neq \text{same reduced J-coupled formula}.
$$

The uncoupled equations stay the same, but the J-coupled factorized formulas must be rederived because the tensor character of $\Omega$ changes the recoupling structure.

If needed, the next step would be to generate AMC output for each of the equations above with `Omega` nonscalar and record the resulting tensor-coupled formulas term by term.

## Explicit AMC probe: exact tensor and scalar intermediates

To fix the phase and index conventions concretely, a focused AMC run was done with input
[amc/examples/sample_input/dc_tensor_intermediate_probe.txt](amc/examples/sample_input/dc_tensor_intermediate_probe.txt)
and output
[amc/examples/sample_output/dc_tensor_intermediate_probe.tex](amc/examples/sample_output/dc_tensor_intermediate_probe.tex).

This probe defines:

- a tensor intermediate from $\Gamma\Omega$: $\Xi^{\beta}$,
- a scalar intermediate from $\Omega\Omega$: $\chi^{\alpha}$,
- a test tensor contraction $f^{probe}$ used only to expose recoupling structure.

AMC returns the following exact reduced forms.

### Tensor intermediate from $\Gamma\Omega$

$$
\Xi^{\beta}_{ij}{}^{\lambda_0}
= \frac12 (-1)^{j_j+\lambda_0}
\sum_{abcJ_0J_1}
(\bar n_a\bar n_b n_c n_i - n_a n_b \bar n_c \bar n_i)
(-1)^{J_0+j_c}
\hat J_0 \hat J_1
\sixj{\lambda_0}{J_1}{J_0}{j_c}{j_i}{j_j}
\Gamma^{J_0J_0 0}_{ciab}
\Omega^{J_0J_1\lambda_0}_{abcj}.
$$

This confirms your rule:

$$
\Gamma\Omega \to \text{tensor intermediate}.
$$

### Scalar intermediate from $\Omega\Omega$

$$
\chi^{\alpha}_{ij}{}^{0}
= \frac12\,\delta_{j_j j_i}\,\hat j_i^{-2}
\sum_{abcJ_0J_1\lambda_0}
(\bar n_a\bar n_b n_c n_i - n_a n_b \bar n_c \bar n_i)
(-1)^{J_0+J_1+\lambda_0}
\hat\lambda_0^{-1}
\Omega^{J_0J_1\lambda_0}_{ciab}
\Omega^{J_1J_0\lambda_0}_{abcj}.
$$

This confirms your second rule in explicit coupled form:

$$
\Omega\Omega \to \text{scalar intermediate}.
$$

### Probe contraction showing tensor recoupling flow

$$
f^{probe}_{ij}{}^{\lambda_0}
= (-1)^{\lambda_0}\hat\lambda_0
\sum_{abJ_0J_1j_0\lambda_1\lambda_2}
(-1)^{j_a+j_b+\lambda_2}
\hat J_0\hat J_1\hat j_0^2
\sixj{j_a}{j_b}{\lambda_1}{j_i}{j_0}{J_0}
\sixj{\lambda_2}{J_0}{J_1}{j_a}{j_j}{j_0}
\sixj{\lambda_0}{\lambda_1}{\lambda_2}{j_0}{j_j}{j_i}
\Xi^{\beta}_{ab}{}^{\lambda_1}
\Omega^{J_0J_1\lambda_2}_{biaj}.
$$

This equation is an explicit recoupling template showing why nonscalar $\Omega$ needs extra angular-momentum sums even
when one keeps the same M-scheme algebra.  It is **not** the stored production $f$; the production scalar equation must
couple/project the relevant tensor pieces back to total rank $0$.
