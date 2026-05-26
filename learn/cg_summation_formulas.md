# Clebsch–Gordan Coefficient Summation Formulas

Reference: **D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii**,
*Quantum Theory of Angular Momentum*, World Scientific, 1988.

---

## Notation

The Clebsch–Gordan (CG) coefficient is written as

$$C^{j_3 m_3}_{j_1 m_1\, j_2 m_2} \;\equiv\; \langle j_1 m_1\, j_2 m_2 \mid j_3 m_3 \rangle$$

It is real in the Condon–Shortley phase convention and vanishes unless

$$m_1 + m_2 = m_3, \qquad |j_1 - j_2| \le j_3 \le j_1 + j_2.$$

The short-hand $\hat{j} \equiv \sqrt{2j+1}$ is used below.

---

## 1. Summation of Two CG Coefficients

### 1.1 Orthogonality I — sum over projections $(m_1, m_2)$

*Varshalovich §8.4, eq. (8.4.3)*

$$\boxed{
\sum_{m_1 m_2}
C^{j m}_{j_1 m_1\, j_2 m_2}\,
C^{j' m'}_{j_1 m_1\, j_2 m_2}
= \delta_{j j'}\,\delta_{m m'}
}$$

**Physical meaning:** The coupled states $\{|j\,m\rangle\}$ form an orthonormal set in $\mathcal{H}_1 \otimes \mathcal{H}_2$.

---

### 1.2 Orthogonality II (Completeness) — sum over $j$ and $m$

*Varshalovich §8.4, eq. (8.4.4)*

$$\boxed{
\sum_{j m}
C^{j m}_{j_1 m_1\, j_2 m_2}\,
C^{j m}_{j_1 m_1'\, j_2 m_2'}
= \delta_{m_1 m_1'}\,\delta_{m_2 m_2'}
}$$

**Physical meaning:** The uncoupled states $\{|j_1 m_1\rangle|j_2 m_2\rangle\}$ are complete in the same tensor-product space (resolution of identity).

---

## 2. Summation of Three CG Coefficients

### 2.1 Sum over $(m_1, m_2)$: recoupling / transformation formula

*Varshalovich §8.7, eq. (8.7.2); also Racah's sum rule*

Summing over $m_1$ and $m_2$ (with $m_{12} = m_1+m_2$, $m = m_{12}+m_3$):

$$\boxed{
\sum_{m_1 m_2}
C^{j_{12}\, m_{12}}_{j_1 m_1\, j_2 m_2}\,
C^{j\, m}_{j_{12}\, m_{12}\, j_3\, m_3}\,
C^{j\, m}_{j_1\, m_1\, j_{23}\, m_{23}}
= (-1)^{j_1+j_2+j_3+j}\,\hat{j}_{12}\,\hat{j}_{23}
\begin{Bmatrix} j_1 & j_2 & j_{12} \\ j_3 & j & j_{23} \end{Bmatrix}
\delta_{m_{23},\,m_2+m_3}
}$$

with $\hat{j}_{12} = \sqrt{2j_{12}+1}$, $\hat{j}_{23} = \sqrt{2j_{23}+1}$, and the braces denote the Wigner 6j symbol.

---

### 2.2 Sum over $(m_1, m_2, m_3)$ producing a 6j symbol

*Varshalovich §8.7, eq. (8.7.3)*

$$\boxed{
\sum_{m_1 m_2 m_3}
(-1)^{j_1-m_1+j_2-m_2+j_3-m_3}\,
C^{j_{12}\, m_{12}}_{j_1 m_1\, j_2 -m_2}\,
C^{j_{23}\, m_{23}}_{j_2 m_2\, j_3 m_3}\,
C^{j_{13}\, m_{13}}_{j_1 m_1\, j_3 -m_3}
= (-1)^{j_1+j_2+j_3}\,\hat{j}_{12}\,\hat{j}_{23}\,\hat{j}_{13}\;
\begin{Bmatrix} j_1 & j_2 & j_{12} \\ j_3 & j_{13} & j_{23} \end{Bmatrix}
(-1)^{j_{12}-m_{12}}\,\delta_{m_{12}+m_{23}-m_{13},\,0}
}$$

---

### 2.3 Symmetric sum over all three projections (Varshalovich §8.7, eq. (8.7.4))

Sum over $m_1, m_2, m_3$ with $m_1+m_2+m_3 = 0$:

$$\boxed{
\sum_{\substack{m_1 m_2 m_3 \\ m_1+m_2+m_3=0}}
C^{0\,0}_{j_1 m_1\, j_2 m_2}\,
C^{j_3 m_3}_{j_1 m_1\, j_2 m_2}\,
C^{j_3 m_3}_{j_1' m_1'\, j_2' m_2'}
= \frac{\delta_{j_3 0}\,\delta_{m_3 0}}{\hat{j}_1}\,\delta_{j_1 j_2}\,\delta_{j_1 j_1'}\,\delta_{j_2 j_2'}
}$$

---

## 3. Supplementary Relations

### 3.1 Relation to Wigner 3j symbols

$$C^{j_3 m_3}_{j_1 m_1\, j_2 m_2}
= (-1)^{j_1-j_2+m_3}\,\hat{j}_3
\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}$$

### 3.2 Racah–Wigner recoupling coefficient

The 6j symbol appears naturally when converting between two different coupling
schemes $(j_1 j_2) j_{12}, j_3 \to J$ and $j_1, (j_2 j_3) j_{23} \to J$:

$$\langle (j_1 j_2) j_{12}, j_3; J \mid j_1, (j_2 j_3) j_{23}; J \rangle
= (-1)^{j_1+j_2+j_3+J}\,\hat{j}_{12}\,\hat{j}_{23}
\begin{Bmatrix} j_1 & j_2 & j_{12} \\ j_3 & J & j_{23} \end{Bmatrix}$$

This is the content of Sections 2.1 and 2.2 above.

---

## 4. Selection Rules and Phase Symmetry

| Property | Formula |
|----------|---------|
| Projection conservation | $C^{jm}_{j_1 m_1\, j_2 m_2} = 0$ unless $m = m_1+m_2$ |
| Triangle rule | $C^{jm}_{j_1 m_1\, j_2 m_2} = 0$ unless $\|j_1-j_2\| \le j \le j_1+j_2$ |
| Exchange symmetry | $C^{jm}_{j_1 m_1\, j_2 m_2} = (-1)^{j_1+j_2-j}\,C^{jm}_{j_2 m_2\, j_1 m_1}$ |
| Time-reversal | $C^{jm}_{j_1 m_1\, j_2 m_2} = (-1)^{j_1+j_2-j}\,C^{j,-m}_{j_1,-m_1\, j_2,-m_2}$ |

---

*See also: Varshalovich et al. (1988), Chapters 8–9; Edmonds (1957), Chapter 3.*
