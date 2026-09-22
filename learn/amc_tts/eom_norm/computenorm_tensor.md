# Tensor `ComputeNorm`: what each diagram is

`ComputeNorm(Qa, Qb)` is \(v_a^{\mathrm{T}} N v_b\).
`ConstructNormMatrix_tensor` builds the bilinear kernel \(N\) for a tensor EOM
operator of rank \(\lambda\) (J), parity \(P\), and isospin rank \(T\).
The numbers that \(N\) multiplies are the **stored** reduced MEs that
`Flatten` / `Unflatten` read (`GetOneBody`, `GetTwoBody` with
`ch_bra ≤ ch_ket`).

The gold check is leftover TTS, not production `Commutator()` on \(T\times T\):

\[
\tfrac12\langle[Q_a,\,Q_b^-]\rangle_\rho
=
\texttt{NormMultiref}(Q_a,Q_b)
=
\tfrac12\,\texttt{GetVSEOM\_Overlap}(\text{leftover }[Q_a,Q_b^-]).
\]

Each diagram is one Wick string of \(\chi^\dagger\rho\chi\) contracted to a
**0-body**. There is no Pandya in \(N\). Angular-momentum 6j / 9j / hats come
from recoupling that string. \(\rho\) is always the **scalar** valence RDM
(`RdmOB`, `RdmTB_J`), even when \(\chi\) is a tensor.

AMC sources and prints: `learn/amc_tts/eom_norm/input/*_tensor.txt` and
`output/*_tensor.tex`. Code: `EOM::ConstructNormMatrix_tensor` in `src/EOM.cc`.
Hats vs leftover overlap: `EQUATION.md`.

---

## How a diagram becomes a number in \(N\)

Four layers. Mixing them is how hats look “special.”

### 1. m-scheme string (physics)

Strip every \(\hat J\) from the scalar kernel. Keep \(\sqrt{2}\) identical-orbit
norms and `Ket::Phase` (those are antisymmetry, not recoupling). Example:
A1 is \(X_{ai}Y_{ai}\), C3 is \(X_{abij}\rho_{abcd}Y_{cdij}\).

### 2. AMC with \(Y\) index-reversed

\(N\) is 0-body (`reduce=false`). \(X\) and \(Y\) are tensors of rank \(\lambda\).
The product \(X_{ai}Y_{ai}\) is not a scalar (M and \(-\mathrm{M}\) do not match).
Reverse the ket so \(\lambda\times\lambda\to 0\): \(X_{ai}Y_{ia}\),
\(X_{abij}Y_{ijab}\). Same pattern as `comm110tts` / `comm220tts`.

AMC then prints a J-scheme formula. Take its **phase, 6j, 9j, and
\(\hat\lambda\) power** as the recoupling. Do **not** take its \(\hat J\)
powers on \(\rho\) as the leftover gold — see § leftover packaging below.

### 3. Fold reversed \(Y\) into the stored ME (always)

IMSRG stores a Hermitian tensor as a reduced RME with `ch_bra ≤ ch_ket`.
The reversed indices AMC wants are the Hermitian partners of what `Flatten`
reads:

| kind | AMC wants | stored | fold in \(N\) |
|---|---|---|---|
| 1b | \(Y_{ia}\) | \(Y_{ai}\) | \(\mathrm{we1}=(-1)^{j_a-j_i}\) |
| 2b | \(Y^{J_1 J_0}_{ijab}\) | \(Y^{J_0 J_1}_{abij}\) | \(\mathrm{we2}=(-1)^{J_1-J_0}\) |

Code: `we1(oa,oi)` and `we2(J0,J1)`. These sit in \(N\) so that
\(v^{\mathrm{T}} N v\) multiplies the same numbers `Flatten` wrote.

### 4. Dagger the **bra** config when the Wick string is \(X^\dagger\)

Same-block diagrams (A1, A2, B2, B4, B5, C1, C3) already put the dagger into
the reversed-\(Y\) fold: the left operator is stored \(X\), the right is stored
\(Y\) after `we1`/`we2`.

Cross-block and pphv leftover strings are written as **daggered bra, stored
ket**. The bra is a 2b config stored \(\langle J_{\mathrm{bra}}\|X\|J_{\mathrm{ket}}\rangle\).
Hermitian dagger:

\[
\langle J_1\|X^\dagger\|J_0\rangle
=
(-1)^{J_1-J_0}\,
\langle J_0\|X\|J_1\rangle.
\]

For a pphv config stored as \(\langle ab\|X\|cd\rangle\) with
\(J_{ab}\) on the particles and \(J_{hv}\) on the hole–valence pair,
dagger \(\to X_{cdab}\) costs

\[
\phi_{\mathrm{dag}}
=
(-1)^{J_{ab}-J_{hv}}.
\]

(`(-1)^{J_{ab}-J_{hv}}=(-1)^{J_{hv}-J_{ab}}\) for integer \(J\).)
This phase is **not** in the AMC print of \(X_{cdab}\). AMC assumes you already
have the daggered ME. The C++ caller multiplies \(\phi_{\mathrm{dag}}\) so it
can keep using the stored \(\langle ab\|X\|cd\rangle\).

1b dagger is the same as `we1`; there is no extra 1b \(\phi_{\mathrm{dag}}\).

---

## Leftover packaging (why C1–C5 hats are not raw AMC)

Gold is leftover commutator \(\times\) `GetVSEOM_Overlap`.
Overlap on a 2b valence pair is \(\rho\cdot\hat J\), not the AMC identity
\(\hat J^2\). Leftover `comm222` / `comm122` already carry
\(\hat J_\rho^{-2}\hat\lambda^{-1}\) (or \(\hat J_{\mathrm{other}}\hat J_\rho^{-1}\hat\lambda^{-1}\)).
So \(N\) must keep **one \(\hat J_\rho\) less** than a naive AMC print:

| leftover | AMC print (typical) | overlap supplies | \(N\) keeps |
|---|---|---|---|
| 110 / 220 identity | \(\hat\lambda^{-1}\) | 0-body m-trace \(\hat\lambda^{-2}\) | \(\hat\lambda^{-2}\) |
| 1b \(\times\rho\) (111 / 221 / 121) | \(\hat\lambda^{-1}\) | \(\rho\cdot\hat\jmath\) in leftover | \(\hat\lambda^{-1}/\hat\jmath\) |
| 2b \(\times\rho\) same block (222, C1/C3) | \(\hat\lambda^{-1}\) | \(\rho\cdot\hat J_\rho\) | \(\hat J_\rho^{-1}\hat\lambda^{-1}\) |
| 2b \(\times\rho\) mixed (122, C2/C5) | \(\hat J_{\mathrm{other}}\hat J_\rho\hat\lambda^{-1}\) | \(\rho\cdot\hat J_\rho\) | \(\hat J_{\mathrm{other}}\hat\lambda^{-1}\) |

These are functions of the **channel \(J\)** and **operator \(\lambda\)** only.
There is no `if (T)` / `if (P)` in the hats. \(T\) and \(P\) only decide which
configs exist (triangle, pair parity, \(|\Delta T_z|=\mathrm{rank}_T\)).

Using \(\hat\lambda^{-2}\) for C1/C3 only matches when \(J_\rho=\lambda\).
That hid the C3 bug at even \(P\), \(T=1\), \(\lambda=3\) (\(J_{pp}=2\)).

---

## Stored configs (so the indices make sense)

| block | `eom_confs` | stored ME |
|---|---|---|
| qv | `{q, v, 0, flat}` | \(X_{qv}\) |
| ph | `{p, h, 0, flat}` | \(X_{ph}\) |
| ppvv | `{ibra, iket, ch_lo, ch_hi}` | \(\langle qq\|X\|vv\rangle\) (ladder drops swapped \(\langle vv\|qq\rangle\)) |
| pphv | same | both \(\langle pp\|X\|hv\rangle\) and \(\langle hv\|X\|pp\rangle\); B1/C4/B3/C5 use **only** \(\langle pp\|hv\rangle\) |
| pphh | same | \(\langle pp\|X\|hh\rangle\) **or** swapped \(\langle hh\|X\|pp\rangle\) (`pphh_legs` maps to logical particles/holes) |

`Jst0`, `Jst1` are the **stored** channel \(J\)s (`ch_lo`, `ch_hi`).
`we2` always uses stored \(J\), not the logical \(J_{pp}\) after a swap.
Identical-orbit \(\sqrt{2}\) and `Ket::Phase` are the scalar-loop factors.

---

## A1 — ph identity (leftover 110)

**In words.** Same particle–hole config on both sides. No \(\rho\). This is
the Fermi piece of \(\tfrac12\langle[Q_a,Q_b^-]\rangle\).

**m-scheme.** \(X_{ai}Y_{ai}\).

**Dagger bra.** No extra \(\phi_{\mathrm{dag}}\). Reversed \(Y\) is the dagger.

**AMC input.** `sum_ai(X_ai * Y_ia)`.

**AMC print.**
\[
N^{(A1)}
=
-\sum_{ai}\,(-1)^{j_a+j_i+\lambda}\,\hat\lambda^{-1}\,
X_{ai}^{\lambda}\,Y_{ia}^{\lambda}.
\]

**Implemented.** Fold \(Y_{ia}=\mathrm{we1}\,Y_{ai}\) and use \(\hat\lambda^{-2}\)
(0-body m-trace):
\[
N_{ii}
=
-(-1)^{j_a+j_i+\lambda}\,\hat\lambda^{-2}\,\mathrm{we1}(a,i).
\]
Same hole on bra and ket (`cf[1]`). Only the diagonal in config space
(identity).

---

## B5 — ph leftover 1b \(\times\rho\) (leftover 111)

**In words.** Same hole, contract the two particles with \(\rho\). Scalar
m-scheme has a leading minus; AMC recoupling cancels it.

**m-scheme.** \(-X_{ai}\rho_{ba}Y_{bi}\).

**Dagger bra.** None beyond `we1`.

**AMC input.** `-sum_abi(X_ai * rho_ba * Y_ib)`.

**AMC print.**
\[
N^{(B5)}
=
\sum_{abi}\,
\delta_{j_b j_a}\,
(-1)^{j_a+j_i+\lambda}\,
\hat\lambda^{-1}\,
X_{ai}^{\lambda}\,\rho_{ba}\,Y_{ib}^{\lambda}.
\]

**Implemented.** Leftover 1b packaging \(\hat\lambda^{-1}/\hat\jmath_a\):
\[
N_{ij}
=
(-1)^{j_a+j_i+\lambda}\,
\hat\lambda^{-1}/\hat\jmath_a\,
\mathrm{we1}(b,i)\,
\rho_{ba},
\]
same hole, \(j_b=j_a\).

---

## B4 — qv leftover 1b \(\times\rho\) (leftover 111) plus 110 identity

**In words.** Same \(q\), contract the two valence legs with \(\rho\).
Also a Fermi 110 on the diagonal (\(n_v(1-n_q)\)).

**m-scheme.** \(X_{ai}\rho_{ij}Y_{aj}\) (same \(q=a\)).

**Dagger bra.** None beyond `we1`.

**AMC input.** `sum_aij(X_ai * rho_ij * Y_ja)`.

**AMC print.**
\[
N^{(B4)}
=
-\sum_{aij}\,
\delta_{j_j j_i}\,
(-1)^{j_a+j_i+\lambda}\,
\hat\lambda^{-1}\,
X_{ai}^{\lambda}\,\rho_{ij}\,Y_{ja}^{\lambda}.
\]

**Implemented.**
\[
N_{ii}^{\mathrm{id}}
=
-(-1)^{j_a+j_i+\lambda}\,\hat\lambda^{-2}\,\mathrm{we1}(a,i)\,n_i(1-n_a),
\]
\[
N_{ij}^{\rho}
=
-(-1)^{j_a+j_i+\lambda}\,
\hat\lambda^{-1}/\hat\jmath_i\,
\mathrm{we1}(a,j)\,
\rho_{ij}.
\]
Same \(q\), \(j_j=j_i\).

---

## A2 — pphh identity (leftover 220)

**In words.** Same pphh config. Fermi \(n_h n_h \bar n_p \bar n_p\).

**m-scheme.** \(X_{abij}Y_{abij}\).

**Dagger bra.** None beyond `we2`. `pphh_legs` still maps swapped
\(\langle hh\|pp\rangle\) to logical particles.

**AMC input.** `sum_abij(X_abij * Y_ijab)`.

**AMC print.**
\[
N^{(A2)}
=
\sum\,
(-1)^{J_0+J_1+\lambda}\,
\hat\lambda^{-1}\,
X^{J_0 J_1\lambda}_{abij}\,
Y^{J_1 J_0\lambda}_{ijab}.
\]

**Implemented.**
\[
N_{ii}
=
(-1)^{J_{\mathrm{st}0}+J_{\mathrm{st}1}+\lambda}\,
\hat\lambda^{-2}\,
\mathrm{we2}(J_{\mathrm{st}0},J_{\mathrm{st}1})\,
n_{h_1}n_{h_2}\bar n_a\bar n_b.
\]

---

## B2 — pphh leftover 1b \(\times\rho\) (leftover 221)

**In words.** Same holes, contract one particle with \(\rho\). Four m-scheme
exchanges of the two particles (not Pandya). Occupied–occupied spectators
have leftover \(N_{\mathrm{occ}}=0\) (He8 valence).

**m-scheme.** \(-X_{abij}\rho_{ca}Y_{cbij}\) plus \(a\leftrightarrow b\),
\(c\leftrightarrow d\).

**Dagger bra.** None beyond `we2`.

**AMC input.** `-sum(X_abij * rho_ca * Y_ijcb)` (one perm; code does all four).

**AMC print.**
\[
N^{(B2)}
=
-\sum\,
\delta_{j_c j_a}\,
(-1)^{J_0+J_1+\lambda}\,
\hat\lambda^{-1}\,
X^{J_0 J_1\lambda}_{abij}\,
\rho_{ca}\,
Y^{J_1 J_0\lambda}_{ijcb}.
\]

**Implemented.** Prefactor
\[
w_0
=
-(-1)^{J_{\mathrm{st}0}+J_{\mathrm{st}1}+\lambda}\,
\hat\lambda^{-1}\,
\mathrm{we2}(J_{\mathrm{st}0},J_{\mathrm{st}1}),
\]
then \(\rho/\hat\jmath\) on the contracted particle, times \(\sqrt{2}\) and
`Ket::Phase` on exchanges. Same holes, same \(J_{pp}\), same \(J_{hh}\).

ppvv has the same four exchanges on the **valence** pair (occupied spectator,
opposite of pphh).

---

## C3 — pphh leftover 2b \(\times\rho\) (leftover 222_pp)

**In words.** Same holes, contract the two particle pairs with \(\rho_2\).
This is the diagram that only fires when both particle pairs are valence and
have even pair parity (two \(p\)-shell orbits). Odd-\(P\) C3 is identically 0
in the He8 \(p\)-shell.

**m-scheme.** \(X_{abij}\rho_{abcd}Y_{cdij}\).

**Dagger bra.** None beyond `we2`.

**AMC input.** `sum(X_abij * rho_abcd * Y_ijcd)`.

**AMC print.**
\[
N^{(C3)}
=
\sum\,
(-1)^{J_0+J_1+J_2}\,
\hat J_2^{-1}\,
X^{J_0 J_1 J_2}_{abij}\,
\rho^{J_0}_{abcd}\,
Y^{J_1 J_0 J_2}_{ijcd}.
\]
(\(J_2=\lambda\).)

**Implemented.** Leftover 222 at \(J_{pp}\) has \(\hat J_{pp}^{-2}\hat\lambda^{-1}\);
overlap supplies \(\rho\cdot\hat J_{pp}\):
\[
N_{ij}
=
(-1)^{J_{\mathrm{st}0}+J_{\mathrm{st}1}+\lambda}\,
\hat J_{pp}^{-1}\hat\lambda^{-1}\,
\mathrm{we2}(J_{\mathrm{st}0},J_{\mathrm{st}1})\,
\rho^{J_{pp}}_{abcd}.
\]
Same holes, \(J_{hh}\), and \(J_{pp}\). Both particle pairs `cvq==1`.
No extra \(T_z\) cut: scalar \(\rho\) already zeros \(\Delta T_z\neq 0\).

---

## C1 — ppvv leftover 2b \(\times\rho\) (leftover 222_pp) plus 220 identity

**In words.** Same \(qq\) (or \(qv\)) spectator, contract the two valence
pairs with \(\rho_2\). Identity only on \(\langle qq\|vv\rangle\), not on
swapped \(\langle vv\|qq\rangle\).

**m-scheme.** \(X_{abij}\rho_{ijkl}Y_{abkl}\).

**Dagger bra.** None beyond `we2`. A \(qv\) spectator (exactly one occupied
leg) flips the leftover-2b sign (`qv_sign = -1`).

**AMC input.** `sum(X_abij * rho_ijkl * Y_klab)`.

**AMC print.**
\[
N^{(C1)}
=
\sum\,
(-1)^{J_0+J_1+J_2}\,
\hat J_2^{-1}\,
X^{J_0 J_1 J_2}_{abij}\,
\rho^{J_1}_{ijkl}\,
Y^{J_1 J_0 J_2}_{klab}.
\]

**Implemented.** Same leftover packaging as C3, \(\rho\) at \(J_{vv}=J_1\):
\[
N_{ij}
=
(\mathrm{qv\_sign})\,
(-1)^{J_0+J_1+\lambda}\,
\hat J_{vv}^{-1}\hat\lambda^{-1}\,
\mathrm{we2}(J_0,J_1)\,
\rho^{J_{vv}}_{ijkl}.
\]
Identity uses \(\hat\lambda^{-2}\) and Fermi \(n_v n_v \bar n_q \bar n_q\).

---

## B1 — pphv leftover 1b \(\times\rho\) (leftover 221) plus 220 identity

**In words.** Same \(pp\) and same hole, contract the two valence legs.
Bra is stored \(\langle ab\|ei\rangle\). The Wick string after dagger is
\(X_{eiab}\,\rho_{ij}\,Y_{abej}\).

**m-scheme.** \(X_{abei}\rho_{ij}Y_{abej}\).

**Dagger bra.** Yes.
\[
\phi_{\mathrm{dag}}=(-1)^{J_{ab}-J_{hv}}.
\]
AMC is written on \(X_{eiab}\) already.

**AMC input.** `sum(X_eiab * rho_ij * Y_abej)`.

**AMC print.**
\[
N^{(B1)}
=
\sum\,
\delta_{j_j j_i}\,
(-1)^{J_0+J_1+\lambda}\,
\hat\lambda^{-1}\,
X^{J_0 J_1\lambda}_{eiab}\,
\rho_{ij}\,
Y^{J_1 J_0\lambda}_{abej}.
\]
Here \(J_0=J_{hv}\) (couples \(e,i\)), \(J_1=J_{ab}\).

**Implemented.** Only \(\langle pp\|hv\rangle\) configs. Leftover 1b packaging:
\[
N_{ij}
=
\phi_{\mathrm{dag}}\,
(-1)^{J_{hv}+J_{ab}+\lambda}\,
\hat\lambda^{-1}/\hat\jmath_c\,
\rho_{c_1 c_2}.
\]
No `we2`: both sides are stored in the same \(\langle pp\|hv\rangle\)
convention, and AMC already used daggered \(X\). Identity is the usual
\(\hat\lambda^{-2}\,\mathrm{we2}\) Fermi factor.

---

## C4 — pphv leftover 2b \(\times\rho\) (leftover 222_ph)

**In words.** Same hole, contract two particles through \(\rho_2\) with the
four scalar permutations of \((e,f)\). This is `Core_Diagram`, not a single
closed AMC of the full C4 string in the kernel loop.

**m-scheme.** \(X_{abcd}\rho_{dfae}Y_{bfce}\) plus three particle exchanges.

**Dagger bra.** Yes. \(\phi_{\mathrm{dag}}=(-1)^{J_{ab}-J_{hv}}\) in the
**caller**, plus `Ket::Phase` on pp swaps (same as the scalar loop).

**AMC of one orbit pair** (`afce_ebfd_tensor`, what `Core_Diagram_tensor`
implements), unreduced scalar \(Z_{abcd}^{J}\):
\[
Z_{abcd}^{J}
=
-\delta_{J J}\,
(-1)^{J+j_a+j_d}
\sum_{ef J_2 J_3 J_4 J_5 j_0}
(-1)^{J_2+J_3+J_4+J_5+j_e+j_f+\lambda}
\hat J_2\hat J_3\hat J_4\hat J_5\,\hat\jmath_0^{2}\,\hat\lambda^{-1}
\]
\[
\times
\begin{Bmatrix} J_3 & \lambda & J_2 \\ j_f & j_a & j_0 \end{Bmatrix}
\begin{Bmatrix} J_4 & \lambda & J_5 \\ j_f & j_d & j_0 \end{Bmatrix}
\begin{Bmatrix}
j_c & j_e & J_3 \\
j_d & J_4 & j_0 \\
J & j_b & j_a
\end{Bmatrix}
X^{J_2 J_3\lambda}_{afce}\,
Y^{J_4 J_5\lambda}_{ebfd}.
\]
Times \(\rho^J\hat J\) and \(\sqrt{2}\), as in scalar `Core_Diagram`.

The full-string AMC `C4_pphv_tensor` (\(X_{cdab}\rho_{dfae}Y_{bfce}\)) is the
same topology with the dagger already in \(X\). The code keeps the scalar
four-perm pipeline and puts \(\phi_{\mathrm{dag}}\) outside
`Core_Diagram_tensor`.

**Implemented.** Only \(\langle pp\|hv\rangle\). Same hole. Four perms of
`(e,f)` with `Ket::Phase`. C4 is the unreduced \(XY\) scalar (\(\hat\lambda^{-1}\)),
not the leftover \(XY-YX\) commutator.

---

## B3 — pphv–ph leftover 1b \(\times\rho\) (leftover 121)

**In words.** Bra is pphv, ket is ph. Always write the **pphv as bra** and
dagger it. Contract the leftover particle with \(\rho\); the hole of the ph
ket is the hole of the pphv.

**m-scheme.** \(X_{cdab}\rho_{bd}Y_{ac}\) (and the other particle:
\(Y_{bc}\,\rho_{ad}\)).

**Dagger bra.** Yes. \(\phi_{\mathrm{dag}}=(-1)^{J_{ab}-J_{hv}}\).

**AMC input.** `sum(X_cdab * rho_bd * Y_ac)`.

**AMC print.**
\[
N^{(B3)}
=
\sum\,
\delta_{j_d j_b}\,
(-1)^{J_{ab}+j_a+j_b}\,
\hat J_{hv}\hat J_{ab}\hat\lambda^{-1}
\begin{Bmatrix} \lambda & J_{ab} & J_{hv} \\ j_b & j_c & j_a \end{Bmatrix}
X^{J_{hv} J_{ab}\lambda}_{cdab}\,
\rho_{bd}\,
Y_{ac}^{\lambda}.
\]
(\(J_0=J_{hv}\), \(J_1=J_{ab}\).)

**Implemented.**
\[
N_{ij}
=
\phi_{\mathrm{dag}}\,
(-1)^{J_{ab}+j_a+j_b}\,
\hat J_{ab}\hat J_{hv}\hat\lambda^{-1}/\hat\jmath_b
\begin{Bmatrix} \lambda & J_{ab} & J_{hv} \\ j_b & j_c & j_a \end{Bmatrix}
\sqrt{2}_{a=b}\,
\rho_{bd}.
\]
The other particle multiplies `Ket::Phase(J_ab)` and \(\rho_{ad}/\hat\jmath_a\).
\(N\) is filled symmetrically (\(N_{ji}=N_{ij}\)).
No `we1` on the ph ket: AMC already used stored \(Y_{ac}\).

---

## C5 — pphv–ph leftover 2b \(\times\rho\) (leftover 122)

**In words.** Same blocks as B3. Contract the two particles of the pphv with
the ph particle and the pphv valence through \(\rho_2\). Scalar m-scheme keeps
the leading minus.

**m-scheme.** \(-X_{cdab}\rho_{abed}Y_{ec}\).

**Dagger bra.** Yes. \(\phi_{\mathrm{dag}}=(-1)^{J_{ab}-J_{hv}}\).

**AMC input.** `-sum(X_cdab * rho_abed * Y_ec)`.

**AMC print.**
\[
N^{(C5)}
=
-\sum\,
(-1)^{J_{ab}+j_d+j_e}\,
\hat J_{hv}\hat J_{ab}\hat\lambda^{-1}
\begin{Bmatrix} J_{hv} & J_{ab} & \lambda \\ j_e & j_c & j_d \end{Bmatrix}
X^{J_{hv} J_{ab}\lambda}_{cdab}\,
\rho^{J_{ab}}_{abed}\,
Y_{ec}^{\lambda}.
\]

**Implemented.** Leftover 122 at \(J_{ab}\) has \(\hat J_{ab}^{-1}\hat J_{hv}\hat\lambda^{-1}\);
overlap supplies \(\rho\cdot\hat J_{ab}\), so drop \(\hat J_{ab}\):
\[
N_{ij}
=
-\phi_{\mathrm{dag}}\,
(-1)^{J_{ab}+j_d+j_e}\,
\hat J_{hv}\hat\lambda^{-1}
\begin{Bmatrix} J_{hv} & J_{ab} & \lambda \\ j_e & j_c & j_d \end{Bmatrix}
\rho^{J_{ab}}_{abed}\,
\sqrt{2}_{e=d}.
\]
Symmetric in \((i,j)\). Both particles of the pphv must be valence.

---

## C2 — ppvv–qv leftover 2b \(\times\rho\) (leftover 122)

**In words.** Bra is ppvv \(\langle ab\|ij\rangle\), ket is qv \(Y_{kb}\).
The spectator particle of the ppvv is the \(q\) of the ket. \(\rho_2\)
contracts the other particle with the valence of the ket, plus the vv pair.

**m-scheme.** \(X_{abij}\rho_{akij}Y_{kb}\).

**Dagger bra.** No 2b \(\phi_{\mathrm{dag}}\) (ppvv is already
\(\langle qq\|vv\rangle\)). Fold the 1b ket: \(Y_{kb}=\mathrm{we1}(b,k)\,Y_{bk}\).

**AMC input.** `sum(X_abij * rho_akij * Y_kb)`.

**AMC print.**
\[
N^{(C2)}
=
\sum\,
(-1)^{J_0+j_a+j_b}\,
\hat J_0\hat J_1\hat\lambda^{-1}
\begin{Bmatrix} J_0 & J_1 & \lambda \\ j_k & j_b & j_a \end{Bmatrix}
X^{J_0 J_1\lambda}_{abij}\,
\rho^{J_1}_{akij}\,
Y_{kb}^{\lambda}.
\]
(\(J_0=J_{qv}\), \(J_1=J_{vv}\).)

**Implemented.** Drop AMC’s extra \(\hat J_{vv}\) (overlap weight):
\[
N_{ij}
=
(-1)^{J_0+j_a+j_b}\,
\hat J_{qv}\hat\lambda^{-1}
\begin{Bmatrix} J_0 & J_1 & \lambda \\ j_k & j_b & j_a \end{Bmatrix}
\mathrm{we1}(b,k)\,
\rho^{J_{vv}}_{akij}\,
\sqrt{2}_{a=k}.
\]
Symmetric in \((i,j)\). The contracted particle \(a\) and both vv legs are
valence.

---

## Phase checklist

| diagram | \(\phi_{\mathrm{dag}}\) on stored bra | fold stored ket | leftover hat in \(N\) |
|---|---|---|---|
| A1 | — | `we1` | \(\hat\lambda^{-2}\) |
| A2 | — | `we2(Jst0,Jst1)` | \(\hat\lambda^{-2}\) |
| B4 id / ρ | — | `we1` | \(\hat\lambda^{-2}\) / \(\hat\lambda^{-1}/\hat\jmath_v\) |
| B5 | — | `we1` | \(\hat\lambda^{-1}/\hat\jmath_p\) |
| B2 | — | `we2` | \(\hat\lambda^{-1}/\hat\jmath_p\) |
| B1 | \((-1)^{J_{ab}-J_{hv}}\) | none (`Y` already stored pp-hv) | \(\hat\lambda^{-1}/\hat\jmath_v\) |
| B3 | \((-1)^{J_{ab}-J_{hv}}\) | none (AMC used \(Y_{ac}\)) | \(\hat J_{ab}\hat J_{hv}\hat\lambda^{-1}/\hat\jmath\) |
| C1 | — | `we2` | \(\hat J_{vv}^{-1}\hat\lambda^{-1}\) |
| C3 | — | `we2` | \(\hat J_{pp}^{-1}\hat\lambda^{-1}\) |
| C2 | — | `we1` on qv | \(\hat J_{qv}\hat\lambda^{-1}\) |
| C5 | \((-1)^{J_{ab}-J_{hv}}\) | none | \(\hat J_{hv}\hat\lambda^{-1}\) |
| C4 | \((-1)^{J_{ab}-J_{hv}}\) in caller | inside `Core_Diagram_tensor` | \(\hat\lambda^{-1}\) on the unreduced \(XY\) |

Universal AMC recoupling phase on same-block 2b leftovers:
\((-1)^{J_0+J_1+\lambda}\). On 1b leftovers:
\((-1)^{j_a+j_i+\lambda}\). That is the \(\lambda\times\lambda\to 0\) couple,
not the Hermitian flip.

---

## What this is not

- \(\rho_3\) is not wired.
- \(\lambda=0\) tensor \(N\) is **not** the scalar `ConstructNormMatrix`
  (reduced vs unreduced hats, and AMC’s extra minus from \(\lambda\times\lambda\to 0\)).
- Gold is leftover TTS (`NormMultiref_tensor`). Production `Commutator()` still
  exits on \(T\times T\).
- Hat powers that differ from AMC are leftover\(\times\)overlap packaging.
  If a mismatch clusters near \(\hat J^{\pm 1}\) or \(1/\sqrt{2J+1}\), it is
  that packaging, not a wrong 6j. See `learn/amc_tts/REDUCED_UNREDUCED.md`.

Checked: random Hermitian pairs, \(v^{\mathrm{T}} N v\) vs leftover
\(\tfrac12\langle[Q_a,Q_b^-]\rangle_\rho\), He8 RDM / He4 Fermi / emax=2,
\(\lambda=1,2,3,4\), \(T=0,1,2\), both parities (24/24).
