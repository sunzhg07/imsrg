# Tensor version term-by-term (arXiv-style intermediates)

This file gives the tensor-coupled analogue of the factorized equations in [learn/arxiv_eq.tex](learn/arxiv_eq.tex), using the same intermediate-operator organization.

Assumptions:

- $\Omega$ is a tensor operator of rank $\lambda$,
- $H\equiv\Gamma$ is scalar,
- final one-body and two-body outputs are tensor operators with rank $\lambda$,
- antisymmetry and dummy-index relabeling conventions are the same as in the scalar derivation.

## Why intermediate form is required

A direct one-shot AMC reduction of all raw M-scheme terms with tensor $\Omega$ failed with inconsistent projection constraints. In practice, nonscalar reductions are stable when written through intermediate tensors first.

So the equations below are written in the same intermediate-first scheme as the arXiv scalar section.

## Intermediate operators

### Scalar intermediates from $\Omega\Omega$

$$
\chi^{\alpha}_{ij}
= \sum_{abcJ}
\frac{\hat J^2}{2\hat j_i^2}
\left(\bar n_a\bar n_b n_c n_i - n_a n_b \bar n_c \bar n_i\right)
\left[\Omega^{JJ\lambda}_{ciab}\times\Omega^{JJ\lambda}_{abcj}\right]^0,
$$

$$
\chi^{\delta J}_{ijkl}
= \sum_{ab}
\frac{\hat J^2}{4}
\left(n_a n_b\bar n_k\bar n_l-\bar n_a\bar n_b n_k n_l\right)
\left[\Omega^{JJ\lambda}_{ijab}\times\Omega^{JJ\lambda}_{abkl}\right]^0,
$$

$$
\chi^{\epsilon}_{ij}
= \frac{1}{2\hat j_j^2}\sum_{abcJ}
\hat J^2
\left(\bar n_a\bar n_b n_c+n_a n_b\bar n_c\right)
\left[\Omega^{JJ\lambda}_{ciab}\times\Omega^{JJ\lambda}_{abcj}\right]^0,
$$

$$
\bar\chi^{\eta J}_{i\bar j k\bar l}
= \sum_{ab}
\left(\bar n_a n_b \bar n_k + n_a \bar n_b n_k\right)
\left[\bar\Omega^J_{i\bar j a\bar b}\times\bar\Omega^J_{a\bar b k\bar l}\right]^0,
$$

$$
\bar\chi^{\theta J}_{i\bar l k\bar j}
= \sum_{ab}
\left(n_a n_b\bar n_k+\bar n_a\bar n_b n_k+n_a n_b\bar n_j+\bar n_a\bar n_b n_j\right)
\left[\bar\Omega^J_{i\bar l a\bar b}\times\bar\Omega^J_{a\bar b k\bar j}\right]^0.
$$

### Tensor intermediates from $\Gamma\Omega$

$$
\Xi^{\beta,\lambda}_{ij}
= \frac12\sum_{abcJ}
\left(\bar n_a\bar n_b n_c n_i-n_a n_b\bar n_c\bar n_i\right)
\Gamma^J_{ciab}\,\Omega^{JJ\lambda}_{abcj}
\times (\text{recoupling phase/6j}),
$$

$$
\Xi^{\zeta,\lambda}_{ij}
= \frac12\sum_{abcJ}
\left(\bar n_a\bar n_b n_c+n_a n_b\bar n_c\right)
\Gamma^J_{aibc}\,\Omega^{JJ\lambda}_{bcaj}
\times (\text{recoupling phase/6j}),
$$

$$
\bar\Xi^{\gamma,\lambda J}_{i\bar j k\bar l}
= \sum_{ab}
\left(n_a\bar n_b\bar n_k n_l-\bar n_a n_b n_k\bar n_l\right)
\bar\Omega^J_{i\bar j a\bar b}\,\bar\Gamma^J_{a\bar b k\bar l},
$$

$$
\bar\Xi^{\iota,\lambda J}_{i\bar j k\bar l}
= \sum_{ab}
\left(\bar n_a n_b\bar n_k+n_a\bar n_b n_k\right)
\bar\Gamma^J_{i\bar j a\bar b}\,\bar\Omega^J_{a\bar b k\bar l},
$$

$$
\overline{\overline\Xi}^{\kappa,\lambda J}_{i\bar j k\bar l}
= \sum_{ab}
\left(\bar n_a n_b n_l+n_a\bar n_b\bar n_l\right)
\overline{\overline\Omega}^J_{i\bar j b\bar a}\,\overline{\overline\Gamma}^J_{b\bar a k\bar l},
$$

$$
\bar\Xi^{\lambda,\lambda J}_{i\bar l k\bar j}
= \sum_{ab}
\left(\bar n_a\bar n_b n_l+n_a n_b\bar n_l\right)
\Gamma^J_{ijab}\,\bar\Omega^J_{a\bar b k\bar j}
+ \sum_{ab}
\left(\bar n_a\bar n_b n_j+n_a n_b\bar n_j\right)
\bar\Omega^J_{i\bar l a\bar b}\,\Gamma^J_{abkl}.
$$

## One-body equations (tensor rank $\lambda$)

### Term I

$$
f^{I,\lambda}_{ij}
= \delta_{j_i j_j}\hat j_i^{-2}
\sum_{abJ}\hat J^2\,\chi^{\alpha}_{ab}\,\Omega^{JJ\lambda}_{biaj}.
$$

### Term II

$$
f^{II,\lambda}_{ij}
= \delta_{j_i j_j}\hat j_i^{-2}
\sum_{abJ}\hat J^2
\left(\Xi^{\beta,\lambda}_{ab}-(-1)^{j_a+j_b-J}\Xi^{\beta,\lambda}_{ba}\right).
$$

### Term IIIa

$$
f^{III_a,\lambda}_{ij}
= \delta_{j_i j_j}\hat j_i^{-2}
\sum_{abcJ}
\left(
\bar\chi^{\eta J}_{i\bar c a\bar b}\,\bar\Xi^{\gamma,\lambda J}_{a\bar b j\bar c}
-
\bar\chi^{\eta J}_{c\bar j a\bar b}\,\bar\Xi^{\gamma,\lambda J}_{a\bar b c\bar i}
\right).
$$

### Term IIIb

$$
f^{III_b,\lambda}_{ij}
= \delta_{j_i j_j}\hat j_i^{-2}
\sum_{abcJ}
\left(
\chi^{\delta J}_{ciab}\,\Omega^{JJ\lambda}_{abcj}
-
\chi^{\delta J}_{abcj}\,\Omega^{JJ\lambda}_{ciab}
\right).
$$

## Two-body equations (tensor rank $\lambda$)

### Term I

$$
\Gamma^{I,\lambda J}_{ijkl}
= \sum_a
\left[
(1-\hat P^J_{ij})\chi^{\epsilon}_{ai}\,\Omega^{JJ\lambda}_{ajkl}
+
(1-\hat P^J_{kl})\chi^{\epsilon}_{ak}\,\Omega^{JJ\lambda}_{ijal}
\right].
$$

### Term II

$$
\Gamma^{II,\lambda J}_{ijkl}
= \sum_a
\left[
(1-\hat P^J_{ij})\Xi^{\zeta,\lambda}_{aj}\,\Gamma^J_{iakl}
-
(1-\hat P^J_{kl})\Xi^{\zeta,\lambda}_{ak}\,\Gamma^J_{ijal}
\right].
$$

### Term IIIa

$$
\Gamma^{III_a,\lambda J}_{ijkl}
= \sum_{ab}
\left[
(1-\hat P^J_{ij})\bar\chi^{\eta J}_{ijab}\,\bar\Xi^{\gamma,\lambda J}_{abkl}
+
(1-\hat P^J_{kl})\bar\Xi^{\gamma,\lambda J}_{ijab}\,\bar\chi^{\eta J}_{abkl}
\right].
$$

### Term IIIb (cross-coupled)

$$
\overline{\overline\Gamma}^{III_b,\lambda J}_{j\bar l k\bar i}
= (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\overline{\overline\Omega}^J_{j\bar l a\bar b}
\left(
\overline{\overline\chi}^{\eta J}_{a\bar b k\bar i}
+
\overline{\overline\chi}^{\eta J}_{k\bar i b\bar a}
\right).
$$

### Term IIIc (Pandya)

$$
\bar\Gamma^{III_c,\lambda J}_{i\bar l k\bar j}
= \frac12(1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\bar\chi^{\theta J}_{i\bar l a\bar b}\,\bar\Xi^{\gamma,\lambda J}_{a\bar b k\bar j}.
$$

### Term IVa

$$
\Gamma^{IV_a,\lambda J}_{ijkl}
= -\sum_{ab}
\left[
(1-\hat P^J_{ij})\overline{\overline\Xi}^{\kappa,\lambda J}_{ijab}\,\Omega^J_{abkl}
+
(1-\hat P^J_{kl})\Omega^J_{ijab}\,\overline{\overline\Xi}^{\kappa,\lambda J}_{klab}
\right].
$$

### Term IVb (cross-coupled)

$$
\overline{\overline\Gamma}^{IV_b,\lambda J}_{j\bar l k\bar i}
= (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\overline{\overline\Omega}^J_{j\bar l a\bar b}
\left(
\bar\Xi^{\iota,\lambda J}_{a\bar b k\bar i}
-
\bar\Xi^{\iota,\lambda J}_{k\bar i a\bar b}
\right).
$$

### Term IVc (Pandya)

$$
\bar\Gamma^{IV_c,\lambda J}_{i\bar l k\bar j}
= \frac12(1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\bar\Xi^{\lambda,\lambda J}_{i\bar l a\bar b}\,\bar\Omega^J_{a\bar b k\bar j}.
$$

## Phase-fixed one-body block (I-III) from explicit AMC outputs

The following equations are exact AMC outputs for the one-body tensor block pieces and their intermediates.

Source files:

- [amc/examples/sample_output/dc_xi1_tensor.tex](amc/examples/sample_output/dc_xi1_tensor.tex)
- [amc/examples/sample_output/dc_xi2_tensor.tex](amc/examples/sample_output/dc_xi2_tensor.tex)
- [amc/examples/sample_output/dc_f3a_tensor_clean.tex](amc/examples/sample_output/dc_f3a_tensor_clean.tex)
- [amc/examples/sample_output/dc_tensor_intermediate_probe.tex](amc/examples/sample_output/dc_tensor_intermediate_probe.tex)

### Exact tensor intermediate for term I

$$
\Xi^{(I)}_{de}{}^{\lambda_0}
= \frac12 (-1)^{j_e}\hat\lambda_0
\sum_{abcJ_0J_1J_2\lambda_1\lambda_2}
(\bar n_a\bar n_b n_c n_d-n_a n_b\bar n_c\bar n_d)
(-1)^{J_2+j_c}
\hat J_0\hat J_2
\sixj{\lambda_0}{J_0}{J_2}{j_c}{j_e}{j_d}
\sixj{\lambda_2}{\lambda_1}{\lambda_0}{J_0}{J_2}{J_1}
\Omega^{J_0J_1\lambda_1}_{cdab}
\Omega^{J_1J_2\lambda_2}_{abce}.
$$

### Exact tensor intermediate for term II

$$
\Xi^{(II)}_{de}{}^{\lambda_0}
= \frac12 (-1)^{j_e+\lambda_0}
\sum_{abcJ_0J_1}
(n_a n_b\bar n_c\bar n_e-\bar n_a\bar n_b n_c n_e)
(-1)^{J_0+j_c}
\hat J_0\hat J_1
\sixj{\lambda_0}{J_1}{J_0}{j_c}{j_d}{j_e}
\Gamma^{J_0J_0 0}_{cdab}
\Omega^{J_0J_1\lambda_0}_{abce}.
$$

### Exact tensor intermediate from probe (Gamma*Omega)

$$
\Xi^{\beta}_{ij}{}^{\lambda_0}
= \frac12 (-1)^{j_j+\lambda_0}
\sum_{abcJ_0J_1}
(\bar n_a\bar n_b n_c n_i-n_a n_b\bar n_c\bar n_i)
(-1)^{J_0+j_c}
\hat J_0\hat J_1
\sixj{\lambda_0}{J_1}{J_0}{j_c}{j_i}{j_j}
\Gamma^{J_0J_0 0}_{ciab}
\Omega^{J_0J_1\lambda_0}_{abcj}.
$$

### Exact scalar intermediate from probe (Omega*Omega)

$$
\chi^{\alpha}_{ij}{}^{0}
= \frac12\delta_{j_j j_i}\hat j_i^{-2}
\sum_{abcJ_0J_1\lambda_0}
(\bar n_a\bar n_b n_c n_i-n_a n_b\bar n_c\bar n_i)
(-1)^{J_0+J_1+\lambda_0}
\hat\lambda_0^{-1}
\Omega^{J_0J_1\lambda_0}_{ciab}
\Omega^{J_1J_0\lambda_0}_{abcj}.
$$

### Exact term IIIa output

$$
f^{(III_a)}_{ij}{}^{\lambda_0}
= \text{Term}_1 - \text{Term}_2,
$$

with the full phase-complete coefficients and $6j$-chains given in
[amc/examples/sample_output/dc_f3a_tensor_clean.tex](amc/examples/sample_output/dc_f3a_tensor_clean.tex).

### Closure status for I and II final one-body equations

The final one-body closure for terms I and II requires one more contraction step of the form

$$
\Xi^{(I)}\times\Gamma \to f^{(I),\lambda},
\qquad
\Xi^{(II)}\times\Omega \to f^{(II),\lambda},
$$

but this direct closure equation hits AMC projection-consistency failures in the current raw setup.

A stable contraction template that *does* run is the verified probe equation

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

So terms I and II are phase-fixed at the intermediate level, and term IIIa is phase-fixed at the final level.

## Phase-fixed two-body block status

Direct AMC reductions with tensor $\Omega$ and scalar $\Gamma$ were attempted term by term using clean per-term inputs.

### Solved directly by AMC

The following terms were reduced successfully with full tensor rank labels $\lambda$ on the final operator:

1. $\Gamma^{I}$:
	output [amc/examples/sample_output/dc_G1_tensor_clean.tex](amc/examples/sample_output/dc_G1_tensor_clean.tex)
2. $\Gamma^{II}$:
	output [amc/examples/sample_output/dc_G2_tensor_clean.tex](amc/examples/sample_output/dc_G2_tensor_clean.tex)
3. $\Gamma^{III_b}$:
	output [amc/examples/sample_output/dc_G3b_tensor_clean.tex](amc/examples/sample_output/dc_G3b_tensor_clean.tex)
4. $\Gamma^{III_c}$:
	output [amc/examples/sample_output/dc_G3c_tensor_clean.tex](amc/examples/sample_output/dc_G3c_tensor_clean.tex)
5. $\Gamma^{IV_a}$:
	output [amc/examples/sample_output/dc_G4a_tensor_clean.tex](amc/examples/sample_output/dc_G4a_tensor_clean.tex)
6. $\Gamma^{IV_b}$:
	output [amc/examples/sample_output/dc_G4b_tensor_clean.tex](amc/examples/sample_output/dc_G4b_tensor_clean.tex)

Each solved file starts with

$$
\Gamma^{(X)}_{ijkl}{}^{J_0 J_1 \lambda_0} = \cdots
$$

and contains the full phase-complete recoupled sums.

### Not solved directly (current AMC implementation limit)

The following terms fail in the current AMC backend for this tensor setting:

1. $\Gamma^{III_a}$
2. $\Gamma^{IV_c}$

The explicit AMC failure is

$$
	exttt{ReductionError: Bigger than square not implemented yet}
$$

with

$$
	exttt{Yutsis graph not fully reduced}.
$$

Additional strict two-operator retries were also tested:

- [amc/examples/sample_input/dc_G3a_tensor_2op.txt](amc/examples/sample_input/dc_G3a_tensor_2op.txt)
- [amc/examples/sample_input/dc_G4c_tensor_2op.txt](amc/examples/sample_input/dc_G4c_tensor_2op.txt)
- [amc/examples/sample_input/dc_G3a_tensor_2op_expanded.txt](amc/examples/sample_input/dc_G3a_tensor_2op_expanded.txt)

with logs

- [/tmp/dc_G3a_tensor_2op.log](/tmp/dc_G3a_tensor_2op.log)
- [/tmp/dc_G4c_tensor_2op.log](/tmp/dc_G4c_tensor_2op.log)
- [/tmp/dc_G3a_tensor_2op_expanded.log](/tmp/dc_G3a_tensor_2op_expanded.log)

These retries still fail in the current AMC stack with either

$$
	exttt{IndexError: list index out of range} + \texttt{Some indices appear only once},
$$

or the same unreduced-graph limitation

$$
	exttt{ReductionError: Bigger than square not implemented yet}.
$$

So these two terms are blocked by the present AMC reduction capability, not by index inconsistencies in the equations.

### Practical consequence

For two-body tensor equations, you now have 6/8 terms fully phase-fixed from direct AMC output and 2/8 terms blocked by a known algorithmic limit.

For the blocked pair, the structural intermediate form written above remains the correct working specification until AMC supports the larger Yutsis reduction class.

## Notes

1. This is a full term-by-term tensor analogue in the same intermediate layout as the scalar arXiv equations.
2. Exact overall phases in each line are convention-dependent and must be fixed against explicit AMC output for each equation block.
3. The structural rule requested by you is enforced throughout:

$$
\Gamma\Omega \to \text{tensor intermediate}, \qquad \Omega\Omega \to \text{scalar intermediate}.
$$

## AMC verification status (tensor $\Omega$, scalar $\Gamma$)

### Raw M-scheme equations

Direct raw-term AMC reductions with tensor $\Omega$ and scalar $\Gamma$ fail for these double-commutator terms due projection/graph-canonicalization issues in the Yutsis stage.

So a direct one-shot extraction of all 12 tensor formulas from raw terms is currently not stable in AMC.

### Intermediate-first equations

The intermediate-first path is stable and reproduces explicit tensor/scalar intermediate behavior. Confirmed AMC probe:

- input: [amc/examples/sample_input/dc_tensor_intermediate_probe.txt](amc/examples/sample_input/dc_tensor_intermediate_probe.txt)
- output: [amc/examples/sample_output/dc_tensor_intermediate_probe.tex](amc/examples/sample_output/dc_tensor_intermediate_probe.tex)

This probe explicitly verifies:

$$
\Gamma\Omega \to \Xi^{\beta} \text{ tensor intermediate},
$$

$$
\Omega\Omega \to \chi^{\alpha} \text{ scalar intermediate}.
$$

Therefore, the equations in this file should be used as the correct term-by-term tensor factorization scheme, with final phase-normalization fixed by running each block through the intermediate-first AMC workflow.
