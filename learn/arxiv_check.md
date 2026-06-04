# AMC check of `arxiv_eq.tex`

## Scope

This note compares the double-commutator equations in [learn/arxiv_eq.tex](learn/arxiv_eq.tex) against explicit angular-momentum reduction with AMC.

The workflow was:

1. take each M-scheme term from the appendix,
2. encode it as an AMC input equation,
3. generate the J-coupled AMC output,
4. compare that output with the J-coupled factorized equations in the arXiv note,
5. allow equivalence under
   - dummy-index relabeling,
   - antisymmetry of two-body matrix elements,
   - Pandya and cross-coupled recoupling.

All AMC runs succeeded for the one-body and two-body pieces.

## Conventions used in the comparison

For antisymmetrized two-body matrix elements in M-scheme,

$$
V_{abcd} = -V_{bacd} = -V_{abdc} = V_{badc}.
$$

In J-scheme, creator-pair and annihilator-pair exchange give the standard phase relation

$$
V^J_{abcd} = -(-1)^{j_a+j_b-J} V^J_{bacd},
$$

$$
V^J_{abcd} = -(-1)^{j_c+j_d-J} V^J_{abdc}.
$$

Dummy summed indices may always be renamed consistently,

$$
\sum_{ab} X_{ab} Y_{ab} = \sum_{xy} X_{xy} Y_{xy},
$$

but slot positions inside barred or double-barred tensors still matter, because they define the coupling scheme.

The Pandya and cross-coupled definitions used in the arXiv note are

$$
\bar A^J_{i \bar l k \bar j} = - \sum_{J'} \hat J'^2
\sixj{j_i}{j_l}{J}{j_k}{j_j}{J'} A^{J'}_{ijkl},
$$

$$
\overline{\overline{A}}^J_{j \bar l k \bar i}
= \sum_{J'} \hat J'^2
\sixj{j_j}{j_l}{J}{j_k}{j_i}{J'}
(-1)^{j_i+j_j-J'} A^{J'}_{ijkl}.
$$

## Files used

The AMC input and output files generated for this check are:

- [amc/examples/sample_input/dc_f1_real.txt](amc/examples/sample_input/dc_f1_real.txt)
- [amc/examples/sample_output/dc_f1_real.tex](amc/examples/sample_output/dc_f1_real.tex)
- [amc/examples/sample_input/dc_f2_real.txt](amc/examples/sample_input/dc_f2_real.txt)
- [amc/examples/sample_output/dc_f2_real.tex](amc/examples/sample_output/dc_f2_real.tex)
- [amc/examples/sample_input/dc_f3a_real.txt](amc/examples/sample_input/dc_f3a_real.txt)
- [amc/examples/sample_output/dc_f3a_real.tex](amc/examples/sample_output/dc_f3a_real.tex)
- [amc/examples/sample_input/dc_f3b_real.txt](amc/examples/sample_input/dc_f3b_real.txt)
- [amc/examples/sample_output/dc_f3b_real.tex](amc/examples/sample_output/dc_f3b_real.tex)
- [amc/examples/sample_input/dc_G1_real.txt](amc/examples/sample_input/dc_G1_real.txt)
- [amc/examples/sample_output/dc_G1_real.tex](amc/examples/sample_output/dc_G1_real.tex)
- [amc/examples/sample_input/dc_G2_real.txt](amc/examples/sample_input/dc_G2_real.txt)
- [amc/examples/sample_output/dc_G2_real.tex](amc/examples/sample_output/dc_G2_real.tex)
- [amc/examples/sample_input/dc_G3a_real.txt](amc/examples/sample_input/dc_G3a_real.txt)
- [amc/examples/sample_output/dc_G3a_real.tex](amc/examples/sample_output/dc_G3a_real.tex)
- [amc/examples/sample_input/dc_G3b_real.txt](amc/examples/sample_input/dc_G3b_real.txt)
- [amc/examples/sample_output/dc_G3b_real.tex](amc/examples/sample_output/dc_G3b_real.tex)
- [amc/examples/sample_input/dc_G3c_real.txt](amc/examples/sample_input/dc_G3c_real.txt)
- [amc/examples/sample_output/dc_G3c_real.tex](amc/examples/sample_output/dc_G3c_real.tex)
- [amc/examples/sample_input/dc_G4a_real.txt](amc/examples/sample_input/dc_G4a_real.txt)
- [amc/examples/sample_output/dc_G4a_real.tex](amc/examples/sample_output/dc_G4a_real.tex)
- [amc/examples/sample_input/dc_G4b_real.txt](amc/examples/sample_input/dc_G4b_real.txt)
- [amc/examples/sample_output/dc_G4b_real.tex](amc/examples/sample_output/dc_G4b_real.tex)
- [amc/examples/sample_input/dc_G4c_real.txt](amc/examples/sample_input/dc_G4c_real.txt)
- [amc/examples/sample_output/dc_G4c_real.tex](amc/examples/sample_output/dc_G4c_real.tex)

## One-body terms

### 1. `f^I_{ij}`

M-scheme term:

$$
f^I_{ij} = \frac12 \sum_{abcde}
(\bar n_a \bar n_b n_c n_d - n_a n_b \bar n_c \bar n_d)
\left( \Omega_{cdab} \Omega_{abce} \Gamma_{eidj}
+ \Omega_{cdab} \Omega_{abce} \Gamma_{diej} \right).
$$

AMC gives the scalar one-body structure

$$
f^I_{ij} = \delta_{j_i j_j} \hat j_i^{-2} \times (\text{J-coupled sum}),
$$

with exactly the two expected operator placements

$$
\Gamma^J_{eidj}, \qquad \Gamma^J_{diej}.
$$

This is consistent with the factorized intermediate

$$
\chi^\alpha_{ij} = \sum_{abcJ}
\frac{\hat J^2}{2\hat j_i^2}
(\bar n_a \bar n_b n_c n_i - n_a n_b \bar n_c \bar n_i)
\Omega^J_{ciab} \Omega^J_{abcj},
$$

and with the arXiv J-coupled factorized form

$$
f^I_{ij} = \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abJ} \hat J^2 \chi^\alpha_{ab} \Gamma^J_{biaj}.
$$

Verdict: consistent.

### 2. `f^{II}_{ij}`

AMC again gives

$$
f^{II}_{ij} = \delta_{j_i j_j} \hat j_i^{-2} \times (\text{J-coupled sum}),
$$

with the two operator placements

$$
\Omega^J_{eidj}, \qquad \Omega^J_{diej}.
$$

The arXiv factorized form

$$
f^{II}_{ij} = \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abJ} \hat J^2 (\chi^\beta_{ab}-\chi^\beta_{ba}) \Omega^J_{biaj}
$$

is the right antisymmetrized J-scheme rewriting. The subtraction

$$
\chi^\beta_{ab}-\chi^\beta_{ba}
$$

is exactly what is needed because swapping the coupled pair in

$$
\Omega^J_{biaj}
$$

produces the antisymmetry phase.

Verdict: consistent.

### 3. `f^{III_a}_{ij}`

AMC produces a recoupled expression with three explicit $6j$ symbols multiplying

$$
\Omega_{abcd} \Omega_{idae} \Gamma_{cejb},
$$

and

$$
\Omega_{abcd} \Omega_{edaj} \Gamma_{cieb}.
$$

This is exactly the pattern that should collapse into the Pandya-coupled factorization

$$
f^{III_a}_{ij} = \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abcJ}
\left(
\bar\chi^{\gamma J}_{i\bar c a\bar b} \Gamma^J_{a\bar b j\bar c}
-
\bar\chi^{\gamma J}_{c\bar j a\bar b} \Gamma^J_{a\bar b c\bar i}
\right).
$$

The barred index placement is consistent with the Pandya definition used in the note.

Verdict: consistent.

### 4. `f^{III_b}_{ij}`

AMC reduces this term without needing Pandya recoupling; the structure matches the standard-coupled intermediate

$$
\chi^{\delta J}_{ijkl} = \sum_{ab}
\frac{\hat J^2}{4}
\left(n_a n_b \bar n_k \bar n_l - \bar n_a \bar n_b n_k n_l \right)
\Omega^J_{ijab} \Omega^J_{abkl}.
$$

The factorized form

$$
f^{III_b}_{ij} = \delta_{j_i j_j} \hat j_i^{-2}
\sum_{abcJ}
\left( \chi^{\delta J}_{ciab} \Gamma^J_{abcj} - \chi^{\delta J}_{abcj} \Gamma^J_{ciab} \right)
$$

matches the index placement seen in the AMC output.

Verdict: consistent.

## Two-body terms

### 5. `\Gamma^I_{ijkl}`

AMC output contains four terms from the two single antisymmetrizers. The operator structure is exactly of the form

$$
(1-\hat P^J_{ij}) \chi^\epsilon_{ai} \Gamma^J_{ajkl}
+ (1-\hat P^J_{kl}) \chi^\epsilon_{ak} \Gamma^J_{ijal}.
$$

This agrees with the arXiv factorized form.

Verdict: consistent.

### 6. `\Gamma^{II}_{ijkl}`

AMC again yields the expected four-term antisymmetrized result, which matches

$$
\Gamma^{II,J}_{ijkl} =
\sum_a \left[
(1-\hat P^J_{ij}) \chi^\zeta_{aj} \Omega^J_{iakl}
-
(1-\hat P^J_{kl}) \chi^\zeta_{ak} \Omega^J_{ijal}
\right].
$$

Verdict: consistent.

### 7. `\Gamma^{III_a}_{ijkl}`

AMC gives the standard-coupled recoupled result expected for the intermediate

$$
\chi^{\eta J}_{ijab}.
$$

The factorized arXiv form

$$
\Gamma^{III_a,J}_{ijkl} =
\sum_{ab}
\left[
(1-\hat P^J_{ij}) \chi^{\eta J}_{ijab} \Gamma^J_{abkl}
+
(1-\hat P^J_{kl}) \Gamma^J_{ijab} \chi^{\eta J}_{abkl}
\right]
$$

is consistent with the AMC index structure after dummy-index relabeling.

Verdict: consistent.

### 8. `\Gamma^{III_b}_{ijkl}`

This term is naturally cross-coupled. AMC expands it into eight antisymmetrized J-coupled contributions. After converting those terms to the double-bar coupling scheme, the natural factorized target is

$$
\overline{\overline{\Gamma}}^{III_b,J}_{j\bar l k\bar i}
=
(1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\overline{\overline\Gamma}^{J}_{j\bar l a\bar b}
\left(
\overline{\overline\chi}^{\eta J}_{a\bar b k\bar i}
+
\overline{\overline\chi}^{\eta J}_{k\bar i b\bar a}
\right).
$$

No contradiction was found between AMC and the arXiv index placement.

Verdict: consistent.

### 9. `\Gamma^{III_c}_{ijkl}`

AMC produces an eight-term antisymmetrized sum. Its recoupled structure is compatible with the single-bar Pandya factorization

$$
\bar\Gamma^{III_c,J}_{i\bar l k\bar j}
=
\frac12 (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\bar\chi^{\theta J}_{i\bar l a\bar b}
\bar\Gamma^J_{a\bar b k\bar j}.
$$

The index placement in the arXiv expression is consistent with the Pandya definition.

Verdict: consistent.

### 10. `\Gamma^{IV_a}_{ijkl}`

AMC gives the correct recoupled J-scheme structure for the M-scheme term

$$
-\sum_{abcd}
(\bar n_c \bar n_d n_a + n_c n_d \bar n_a)
\left[(1-\hat P_{ij}) \Omega_{aicd} \Omega_{dbkl} \Gamma_{jcba}
+
(1-\hat P_{kl}) \Omega_{dcak} \Omega_{ijcb} \Gamma_{bald}
\right].
$$

The J-coupled arXiv line is

$$
\Gamma^{IV_a,J}_{ijkl} =
-\sum_{ab}
\left[
(1-\hat P^J_{ij}) \chi^{\kappa J}_{ijab} \Omega^J_{abkl}
+
(1-\hat P^J_{kl}) \Omega^J_{ijab} \chi^{\kappa J}_{klab}
\right].
$$

However, the earlier uncoupled factorized line in the note uses

$$
\Gamma^{IV_a}_{ijkl} =
-\sum_{ab}
\left[
(1-\hat P_{ij}) \chi^\kappa_{ijab} \Omega_{bakl}
-
(1-\hat P_{kl}) \Omega_{ijab} \chi^\kappa_{klba}
\right].
$$

<span style="color:red">Suspicious:</span> the uncoupled factorized form and the J-coupled factorized form do not use the same `ab` ordering.

The mismatch is

$$
\Omega_{bakl} \leftrightarrow \Omega^J_{abkl},
$$

$$
\chi^\kappa_{klba} \leftrightarrow \chi^{\kappa J}_{klab}.
$$

This may still be repairable by antisymmetry phases and consistent relabeling, but it is not manifestly identical as written.

Verdict: <span style="color:red">needs explicit cleanup in the manuscript</span>.

### 11. `\Gamma^{IV_b}_{ijkl}`

AMC gives an eight-term cross-coupled result. The arXiv factorized target,

$$
\overline{\overline\Gamma}^{IV_b,J}_{j\bar l k\bar i}
=
(1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\overline{\overline\Omega}^{J}_{j\bar l a\bar b}
\left(
\bar\chi^{\iota J}_{a\bar b k\bar i}
-
\bar\chi^{\iota J}_{k\bar i a\bar b}
\right),
$$

is compatible with the AMC operator slots after cross-coupling.

Verdict: consistent.

### 12. `\Gamma^{IV_c}_{ijkl}`

AMC gives an eight-term Pandya-type result consistent with the M-scheme term and with a barred factorization.

The arXiv expression is

$$
\bar\Gamma^{IV_c,J}_{i\bar l k\bar j}
=
\frac12 (1-\hat P^J_{ij})(1-\hat P^J_{kl})
\sum_{ab}
\bar\chi^{\lambda J}_{i\bar l a\bar b}
\bar\Omega^J_{a\bar k j\bar b}.
$$

<span style="color:red">Suspicious:</span> the barred Omega index placement is not the most natural one from the Pandya definition and from the uncoupled factorized operator order.

The current written factor is

$$
\bar\Omega^J_{a\bar k j\bar b}.
$$

A more natural Pandya target from the operator ordering is of the form

$$
\bar\Omega^J_{a\bar b k\bar j}
\quad \text{or an equivalent relabeled version in the same coupling scheme,}
$$

not one that mixes the external `j` into the third slot while the internal `b` sits in the fourth slot without an explicit derivation.

Verdict: <span style="color:red">index placement should be rederived carefully</span>.

## Intermediate checks

The one-body intermediates

$$
\chi^\alpha, \chi^\beta, \bar\chi^{\gamma J}, \chi^{\delta J}
$$

are consistent with the AMC outputs.

The two-body intermediates

$$
\chi^\epsilon, \chi^\zeta, \chi^{\eta J}, \bar\chi^{\theta J}, \bar\chi^{\iota J}
$$

are also consistent with the AMC reductions.

For the remaining two-body intermediates, the main concerns are presentational consistency rather than an AMC failure:

<span style="color:red">Suspicious:</span> the J-coupled definition of `\chi^\lambda` is missing the superscript `J` in the note.

Current text:

$$
\chi^{\lambda~}_{ijkl} = \cdots
$$

Expected notation:

$$
\chi^{\lambda J}_{ijkl} = \cdots
$$

Also, the note defines

$$
\overline{\overline\chi}^{\kappa J}_{i\bar j k\bar l}
$$

while the factorized `IV_a` line uses the standard-coupled object

$$
\chi^{\kappa J}_{ijab}.
$$

This is not automatically wrong, but it means an implicit recoupling step is being used and should be made explicit.

## Pandya and cross-coupling summary

The arXiv note is using the correct general strategy:

- `III_a` and `III_c` are naturally expressed with single-bar Pandya objects,
- `III_b` and `IV_b` are naturally expressed with double-bar cross-coupled objects,
- `I`, `II`, `III_a`, `IV_a` can be written in standard coupling after defining suitable intermediates.

The AMC outputs confirm that this coupling-scheme split is physically and algebraically sensible.

The remaining issues are not with AMC. They are local index-order and notation issues inside the manuscript's factorized J-coupled presentation.

## Bottom line

The M-scheme equations reduce successfully with AMC, and most of the J-coupled factorized equations in the note are consistent with the AMC reductions.

The places that should be checked and probably corrected in the arXiv note are:

- <span style="color:red">`\Gamma^{IV_a}`: mismatch between `ab` and `ba` ordering across the uncoupled and J-coupled factorized forms.</span>
- <span style="color:red">`\Gamma^{IV_c}`: suspicious barred Omega index order in the J-coupled factorized form.</span>
- <span style="color:red">`\chi^{\lambda}`: missing J label in the J-coupled intermediate definition.</span>
- <span style="color:red">`\chi^{\kappa}`: coupling-scheme translation is implicit and should be stated explicitly.</span>

Everything else checked here is consistent with AMC after allowing dummy-index relabeling, antisymmetry phases, and Pandya/cross-coupled recoupling.
