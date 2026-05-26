# `EOM::ArnoldiSolve` — algorithm explained

This note explains the algorithm in
[src/EOM.cc](../src/EOM.cc#L2076).

We solve an eigenvalue problem for
$$ H \;=\; H_1 + H_2 $$
where:

- $H_1 v$ is **affordable** (cheap one-sided action).
- $H_2 v$ is **not affordable**; only diagonal expectation values
  $\langle u|H_2|u\rangle$ are.
- The full $H$ is **symmetric** in the metric
  $\langle a|b\rangle = \texttt{ComputeNorm}(a,b)$, but $H_1$ and $H_2$
  individually are **not** symmetric.

---

## TL;DR — equivalence to the original code

The new code computes **the exact same numbers** as the original. The
only difference is the diagonal $i = j$ short-circuit, which is an
algebraic identity (proved in §5), not a new formula.

| Quantity | Original code | New code |
|---|---|---|
| Off-diagonal $H_{ij}$ ($i \neq j$) | $s_1 + s_2$ | $s_1 + s_2$ (identical) |
| Diagonal $H_{jj}$ | $\tfrac12(4 - 2)\,h_2^\text{diag}[j]$ via polarization on $2 v_j$ | $h_2^\text{diag}[j]$ directly from cache |
| Eigensolve | `eig_sym(hall_sub)` | `eig_sym(symmatu(0.5*(sub+sub^T)))` |

The diagonal change saves one expensive `DcomMultiref` per step and
avoids a $4x - 2x \to x$ floating-point cancellation. The eigensolve
change forces the FP-noise asymmetry of `hall` to zero before
diagonalization. Neither changes any matrix element.

Everything below derives why $H_{ij} = s_1 + s_2$ is correct from the
symmetry assumption.

---

## 1. The symmetry assumption (your input)

You told the code that $H = H_1 + H_2$ is symmetric:
$$
\langle i|H_1|j\rangle + \langle i|H_2|j\rangle
\;=\;
\langle j|H_1|i\rangle + \langle j|H_2|i\rangle.
\tag{S}
$$

For each operator $A \in \{H_1, H_2\}$, split its matrix element into
symmetric and antisymmetric pieces:
$$
\langle i|A|j\rangle = s_A + a_A, \qquad
\langle j|A|i\rangle = s_A - a_A,
$$
with
$$
s_A = \tfrac12\bigl(\langle i|A|j\rangle + \langle j|A|i\rangle\bigr),
\qquad
a_A = \tfrac12\bigl(\langle i|A|j\rangle - \langle j|A|i\rangle\bigr).
$$

Plugging into (S):
$$
(s_1 + a_1) + (s_2 + a_2) \;=\; (s_1 - a_1) + (s_2 - a_2)
\quad\Longrightarrow\quad
\boxed{\,a_1 + a_2 = 0\,}.
\tag{1}
$$

So the antisymmetric parts of $H_1$ and $H_2$ are exact opposites.
Therefore
$$
H_{ij}
= \langle i|H_1|j\rangle + \langle i|H_2|j\rangle
= (s_1 + a_1) + (s_2 + a_2)
\stackrel{(1)}{=} s_1 + s_2.
\tag{2}
$$

**This is the central identity used throughout.**

---

## 2. Computing $s_1$ — the H1 piece

We have one-sided $H_1 v$ (cached as `h1v_cache[k]`). Two inner
products give us
$$
\langle i|H_1|j\rangle = \langle v_i \,|\, H_1 v_j\rangle, \qquad
\langle j|H_1|i\rangle = \langle v_j \,|\, H_1 v_i\rangle,
$$
and we symmetrize:
$$
\boxed{\,s_1 \;=\; \tfrac12\bigl(\langle v_i|H_1 v_j\rangle + \langle v_j|H_1 v_i\rangle\bigr)\,}.
$$

In the code:
```cpp
double h1ij = ComputeNorm(lanczos_vector[i], h1v_j);
double h1ji = ComputeNorm(lanczos_vector[j], h1v_cache[i]);
h1_sym      = 0.5 * (h1ij + h1ji);
```

---

## 3. Computing $s_2$ — the H2 piece via polarization

Expand the cross term as you wrote:
$$
\langle u+v|H_2|u+v\rangle
= \langle u|H_2|u\rangle + \langle v|H_2|v\rangle
+ \langle u|H_2|v\rangle + \langle v|H_2|u\rangle.
$$

Rearrange:
$$
\tfrac12\bigl(\langle u+v|H_2|u+v\rangle - \langle u|H_2|u\rangle - \langle v|H_2|v\rangle\bigr)
= \tfrac12\bigl(\langle u|H_2|v\rangle + \langle v|H_2|u\rangle\bigr)
= s_2.
$$

This polarization gives us **only the symmetric part $s_2$** — the
antisymmetric part $a_2$ cannot be extracted from diagonal
expectation values alone, exactly as you noted. **But that is fine**,
because by identity (2) the matrix element we actually want is
$s_1 + s_2$, and $s_2$ is precisely what polarization delivers.

So
$$
\boxed{\,s_2 \;=\; \tfrac12\bigl(h_2^\text{cross}[i,j] - h_2^\text{diag}[i] - h_2^\text{diag}[j]\bigr)\,},
$$
with $h_2^\text{cross}[i,j] = \langle v_i + v_j|H_2|v_i + v_j\rangle$.

In the code:
```cpp
Operator v_sum  = lanczos_vector[i] + lanczos_vector[j];
double h2_cross = DcomMultiref(Hs, v_sum).first;
h2_sym          = 0.5 * (h2_cross - h2_diag[i] - h2_diag[j]);
```

We cache $h_2^\text{diag}[k] = \langle v_k|H_2|v_k\rangle$ once per
basis vector (as soon as $v_k$ is added), so each off-diagonal needs
only **one new cross-term evaluation**.

---

## 4. Putting it together — the off-diagonal formula

By (2), §2, and §3, for every $i \neq j$:
$$
\boxed{\,H_{ij} \;=\; s_1 + s_2\,}, \qquad H_{ji} = H_{ij}.
$$

Code:
```cpp
hall(i, j) = hall(j, i) = h1_sym + h2_sym;
```

This is **identical** to the original implementation.

---

## 5. The diagonal short-circuit ($i = j$) — algebraic identity

What happens if we naively reuse the polarization formula at $i = j$?
Set $u = v = v_j$, so $u + v = 2 v_j$. By bilinearity of $H_2$,
$$
\langle 2v_j|H_2|2v_j\rangle = 4\,\langle v_j|H_2|v_j\rangle = 4\,h_2^\text{diag}[j],
$$
hence
$$
\underbrace{\tfrac12\bigl(\langle 2v_j|H_2|2v_j\rangle - h_2^\text{diag}[j] - h_2^\text{diag}[j]\bigr)}_{\text{old code}}
= \tfrac12(4 - 2)\,h_2^\text{diag}[j]
= \underbrace{h_2^\text{diag}[j]}_{\text{new code}}.
$$

The two are **algebraically identical**. The new code uses the
right-hand side directly:

```cpp
if (i == j) {
    h1_sym = ComputeNorm(lanczos_vector[j], h1v_j);  // <v|H1|v>
    h2_sym = h2_diag[j];                              // <v|H2|v>
}
```

Pure optimization:

- **Cost:** saves one `DcomMultiref(Hs, 2 v_j)` per step (the most
  expensive call in the inner loop).
- **Numerics:** avoids the cancellation $\tfrac12(4x - 2x) = x$, which
  loses ~1 decimal digit when $x = h_2^\text{diag}[j]$ is large.
- **Result:** identical to the original up to FP noise the original
  itself suffered from.

---

## 6. Why the antisymmetric parts must really cancel

Suppose your $H_1/H_2$ split is *almost* but not exactly consistent
with (S), i.e. $a_1 + a_2 = \varepsilon \neq 0$. Then:

- The polarization on $u+v$ still returns exactly $s_2$ (it cannot
  return anything else; $a_2$ is invisible to it by construction).
- The H1 symmetrization still returns exactly $s_1$.
- So `hall(i,j)` is still $s_1 + s_2$.

But the **true** $\langle i|H|j\rangle = s_1 + s_2 + \varepsilon$
while $\langle j|H|i\rangle = s_1 + s_2 - \varepsilon$. The two are
not equal, $H$ is genuinely asymmetric, and `eig_sym(hall)` is
solving the wrong problem. The code prints a diagnostic so you can
detect this case:
$$
r_{ij} \;=\; \frac{|\langle v_i|H_1 v_j\rangle - \langle v_j|H_1 v_i\rangle|}{2\,|H_{ij}|}
       \;=\; \frac{|a_1|}{|s_1 + s_2|}.
$$
A warning fires when $\max_{i<j} r_{ij} > 10^{-3}$. This is
informational — large $|a_1|$ is allowed if $|a_2| = |a_1|$ matches
it, but it tells you whether your split is internally consistent.

---

## 7. Subspace generation (unchanged)

The Lanczos basis is built from $H_1$ alone:

1. Normalize $v_0 = v_i / \sqrt{\langle v_i|v_i\rangle}$.
2. At step $j$, form $w \;=\; H_1 v_j$ (`HtcMultiref(Hs, v_j)`).
3. Apply `ProjectOprator(w)` to enforce subspace constraints.
4. Double-pass classical Gram–Schmidt against
   $\{v_0, \dots, v_j\}$ in the metric `ComputeNorm`, with another
   `ProjectOprator` between passes.
5. $v_{j+1} \;=\; w / \sqrt{\langle w|w\rangle}$.

The resulting subspace is the Krylov space $\mathcal K_m(H_1, v_0)$,
**not** $\mathcal K_m(H, v_0)$. We never apply $H_2$ to a vector.

---

## 8. Diagonalization with explicit symmetrization

Even though we write `hall(i,j) = hall(j,i) = value`, repeated
read–modify–write on the matrix can leave bit-level FP asymmetry that
confuses `arma::eig_sym`. Every diagonalization is routed through

```cpp
auto eigensolve_sub = [&hall](int dim, arma::vec &ev, arma::mat &evec) {
    arma::mat sub = hall.submat(0, 0, dim - 1, dim - 1);
    sub = arma::symmatu(0.5 * (sub + sub.t()));
    arma::eig_sym(ev, evec, sub);
};
```

In exact arithmetic this is a no-op; in FP it removes a class of
subtle bugs.

---

## 9. Convergence and breakdown handling (unchanged logic, factored code)

After each new basis vector, if $j+1 \ge \texttt{min\_iter}$:

- Diagonalize the leading $(j{+}1)\times(j{+}1)$ block of `hall`.
- Take the lowest `state_want` eigenvalues as the current Ritz
  estimates `e`.
- If $\max|e - \texttt{prev\_e}| < \texttt{tol}$ (absolute or
  relative), declare convergence and stop.

Three breakdown modes can stop the iteration early; each one solves
the current subspace via `eigensolve_sub` and returns:

| Condition | Meaning |
|---|---|
| $\lvert b_j \rvert < \texttt{null\_tol} \cdot \texttt{cn0}$ | New vector exhausted by projection (null space). |
| $\lvert b_j \rvert < \texttt{bj\_tol}$ | Exact breakdown: $w \approx 0$. |
| $b_j < 0$ | Indefinite metric: `ComputeNorm` returned a negative norm. |

A final eigensolve over the entire accumulated subspace runs
unconditionally before returning.

---

## 10. Variational consistency check (new diagnostic only)

Every 5 steps the code now also computes, for each Ritz vector
$|\Psi_k\rangle = \sum_m c_{mk}\,v_m$,
$$
\mathcal E(\Psi_k) \;=\; \frac{\langle\Psi_k|H_1|\Psi_k\rangle + \langle\Psi_k|H_2|\Psi_k\rangle}{\langle\Psi_k|\Psi_k\rangle},
$$
via `EOM::ExpectationValue(Psi_k)`, and prints it next to the Ritz
value $E_k$ from `hall`.

In **exact arithmetic** the two must agree, provided
$\langle v_i|v_j\rangle = \delta_{ij}$ and `hall(i,j)` is correctly
$\langle v_i|H|v_j\rangle$. A growing discrepancy indicates either
loss of orthogonality or a linearity / projection issue inside
`HtcMultiref` / `DcomMultiref` when called on a generic superposition.

This is purely a diagnostic — it does not change the algorithm.

---

## 11. Ritz vectors

For each $E_k$ with eigenvector $\mathbf c_k = (c_{0k}, c_{1k}, \dots)^T$:
$$
|\Psi_k\rangle \;=\; \sum_{m=0}^{n_b-1} c_{mk}\, v_m,
$$
assembled in the loop that fills `ritz_vecs`.

---

## 12. Cost summary per step $j$

| Quantity | Count | Cost class |
|---|---|---|
| $H_1 v_j$ | 1 | one-sided H1 action |
| $\langle v_i \mid H_1 v_j\rangle, \langle v_j \mid H_1 v_i\rangle$ | $2j$ inner products | cheap |
| $\langle v_i+v_j \mid H_2 \mid v_i+v_j\rangle$, $i < j$ | $j$ | one diagonal H2 expectation each |
| $\langle v_{j+1} \mid H_2 \mid v_{j+1}\rangle$ (cache) | 1 | one diagonal H2 expectation |
| Diagonal $\langle v_j \mid H_2 \mid v_j\rangle$ | **0** (uses cache; old code did 1) | — |
| GS reorthogonalization | $2(j+1)$ inner products + projections | cheap |
| Eigensolve of $(j{+}1)\times(j{+}1)$ | 1 | tiny dense |

Total expensive H2 evaluations after $m$ steps:
$\binom{m}{2} + m = m(m+1)/2 = \mathcal O(m^2)$,
**one fewer per step than the original** (which also paid for the
diagonal polarization), with no $H_2 v$ ever required.
