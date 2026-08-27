# AMC from BruteForce-mapped m-scheme vs scalar direct

**Inputs** (P stripped; code restores \(P_{pg}P_{qh}\)):

| File | M-scheme |
|---|---|
| `input/IIb_mscheme.txt` | \(-\sum(\bar n_b n_c n_d+\bar n_c\bar n_d n_b)\,\Omega_{dcbk}\Omega_{biac}\Gamma_{jald}\) |
| `input/IId_mscheme.txt` | \(-\sum(\bar n_c n_b n_d+\bar n_b\bar n_d n_c)\,\Omega_{jcbd}\Omega_{balc}\Gamma_{diak}\) |
| `input/G3b_IIb_IId_mscheme.txt` | IIb + IId |

**Outputs:** `output/{IIb,IId,G3b_IIb_IId}_mscheme.tex`

## Verdict: **MATCH** scalar BruteForce IIb / IId

Index map \(i,j,k,l \leftrightarrow p,g,q,h\), \(\delta_{J_0 J_1}\).

### IIb

| Piece | AMC | BruteForce L17307–17312 / L17412 |
|---|---|---|
| occ | \(\bar n_b n_c n_d+\bar n_c\bar n_d n_b\) | same |
| hats | \(\hat J_2^2\hat J_3^2\hat J_4^2\hat J_5^2\) | \((2J+1)\) four times |
| 6j | \(\{j_k j_d J_5;\,j_c j_b J_2\}\cdots\{J_0 J_4 J_5;\,j_a j_i j_j\}\) | \(\{j_q j_d J_5;\,j_c j_b J_2\}\cdots\{J_0 J_4 J_5;\,j_a j_p j_g\}\) |
| ops | \(\Omega_{dcbk}\Omega_{biac}\Gamma_{jald}\) | \(\eta_{dcbq}\eta_{bpac}\Gamma_{gahd}\) |
| sign | \(-\) | \(-=\) |

### IId

| Piece | AMC | BruteForce L17552–17556 / L17659 |
|---|---|---|
| occ | \(\bar n_c n_b n_d+\bar n_b\bar n_d n_c\) | same (code \(+\); comment typo \(-\)) |
| 6j | \(\{j_d j_j J_5;\,\ldots\}\{J_5 J_4 J_0;\,j_k j_l j_a\}\) | \(\{j_d j_g J_5;\,\ldots\}\{J_5 J_4 J_0;\,j_q j_h j_a\}\) |
| ops | \(\Omega_{jcbd}\Omega_{balc}\Gamma_{diak}\) | \(\eta_{gcbd}\eta_{bahc}\Gamma_{dpaq}\) |

So the m-scheme mapped from the scalar direct code **does** regenerate the scalar JT equations in `ReferenceImplementations.cc`.
