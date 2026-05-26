# Arxiv LaTeX equations

This file reproduces the LaTeX equation blocks from `learn/arxiv.txt` as display math suitable for Markdown renderers with MathJax/KaTeX support.



The one-body piece of the double commutator is
$$
\begin{align}
%%%%%%%%%%%%%%%%%  I
 f^I_{ij}  =&\frac{1}{2} \sum_{abcde} \left( \bar{n}_a \bar{n}_b  {n}_c  {n}_d  -  {n}_a  {n}_b \bar{n}_c  \bar{n}_d  \right)   \nonumber \\
  &\times \left( \Omega_{cdab}\Omega_{abce}   \Gamma_{eidj}  +  \Omega_{cdab} \Omega_{abce}  \Gamma_{diej}  \right)\\
 %%%%%%%%%%%%%%%%  II   
  f^{II}_{ij}  =& \frac{1}{2} \sum_{abcde} \left( {n}_a  {n}_b  \bar{n}_c  \bar{n}_e  -   \bar{n}_a  \bar{n}_b  {n}_c  {n}_e  \right) \nonumber \\
 &\times  \left( \Gamma_{cdab} \Omega_{abce}   \Omega_{eidj}   - \Gamma_{cdab} \Omega_{abce}  \Omega_{diej}  \right)  \\
 %
f^{III_a}_{ij} =& \sum_{abcde} \left( \bar{n}_a  \bar{n}_b  {n}_c  {n}_d  -  {n}_a  {n}_b  \bar{n}_c  \bar{n}_d \right) \nonumber \\
 &\times   \left( \Omega_{abcd}  \Omega_{idae} \Gamma_{cejb}  -  \Omega_{abcd}  \Omega_{edaj} \Gamma_{cieb}  \right)\\
  %
f^{III_b}_{ij}  =& \frac{1}{4} \sum_{abcde} \left( \bar{n}_a  \bar{n}_b  {n}_c {n}_d  -  {n}_a  {n}_b  \bar{n}_c  \bar{n}_d \right) \nonumber \\
 &\times   \left( \Omega_{abcd} \Omega_{cdej} \Gamma_{eiab} -  \Omega_{abcd} \Omega_{eiab} \Gamma_{cdej} \right)
    %
\end{align}
$$

These may be re-expressed in a factorized form
$$
\begin{align}  
 f^I_{ij}  =& \sum_{ab} \left( \chi^{\alpha}_{ab}  \Gamma_{biaj} + \chi^{\alpha}_{ab} \Gamma_{aibj} \right) \\
   %%%%%%%%%%
f^{II}_{ij}   =& \sum_{ab} \left( \chi^{\beta}_{ab} \Omega_{biaj} - \chi^{\beta}_{ab} \Omega_{aibj} \right)\\
  %%%%%%%%%%
 f^{III_a}_{ij}  =& \sum_{abc} \left( \chi^{\gamma}_{icab} \Gamma_{abjc} -  \Gamma_{ciab}  \chi^{\gamma}_{abcj} \right)\\
  %%%%%%%%%%
f^{III_b}_{ij}  =& \sum_{abc} \left(\chi^{\delta}_{ciab}  \Gamma_{abcj} - \Gamma_{ciab} \chi^{\delta}_{abcj} \right)
\end{align}
$$
where we have defined the intermediate quantities
$$
\begin{align} 
\chi^{\alpha}_{ij} =&\frac{1}{2} \sum_{abc} \left( \bar{n}_a \bar{n}_b  {n}_c  {n}_i  -  {n}_a {n}_b \bar{n}_c  \bar{n}_i \right)  \Omega_{ciab} \Omega_{abcj}  \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chi^{\beta}_{ij} =&\frac{1}{2} \sum_{abc} \left(  \bar{n}_a   \bar{n}_b {n}_c {n}_i  -  {n}_a  {n}_b   \bar{n}_c   \bar{n}_i \right) \Omega_{ciab} \Gamma_{abcj} \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chi^{\gamma}_{ijkl} =& \sum_{ab} \left(  {n}_a \bar{n}_b {n}_j  \bar{n}_k   -  \bar{n}_a {n}_b \bar{n}_j  {n}_k   \right)  \Omega_{ajkb} \Omega_{ibal} \\
% 
\chi^{\delta}_{ijkl} =& \frac{1}{4}\sum_{ab} \left( \bar{n}_i \bar{n}_j {n}_a {n}_b  -  {n}_i  {n}_j  \bar{n}_a \bar{n}_b  \right)  \Omega_{ijab} \Omega_{abkl}
\end{align}
$$
We use superscript Greek letters to distinguish the intermediates.

The two-body piece of double the commutator is
$$
\begin{align}
%  Ia
\Gamma^{I}_{ijkl}  =&\frac{1}{2} \sum_{abcd} \left( \bar{n}_a  \bar{n}_b  {n}_c  +  {n}_a  {n}_b  \bar{n}_c \right) \nonumber\\
& \times \left\{ (1\! -\! \hat P_{ij} ) \Omega_{ciab} \Omega_{abcd} \Gamma_{djkl} \right. \nonumber\\
&~+ \left.  (1\! -\! \hat P_{kl} ) \Omega_{cdab} \Omega_{abcl} \Gamma_{ijkd} \right\}\\
%    IV
\Gamma^{II}_{ijkl}   =& 
- \frac{1}{2} \sum_{abcd} \left(   \bar{n}_a \bar{n}_b {n}_c +  \bar{n}_c  {n}_a  {n}_b \right)  \nonumber \\
&\times  \left\{  (1\! -\! \hat P_{ij} )  \Omega_{cjab}  \Gamma_{abcd}  \Omega_{idkl}\right. \nonumber\\
&~+ \left.  (1\! -\! \hat P_{kl} ) \Gamma_{cdab}  \Omega_{abcl} \Omega_{ijkd}  \right\}\\
%   IIa
\Gamma^{III_a}_{ijkl} =& - \sum_{abcd} \left( \bar{n}_c  \bar{n}_d  {n}_a  +  {n}_a  \bar{n}_c  \bar{n}_d \right)  \nonumber\\
& \times \left\{  (1\! -\! \hat P_{ij} )  \Omega_{ajcd} \Omega_{idab} \Gamma_{cbkl} \right. \nonumber\\
& ~+ \left. (1\! -\! \hat P_{kl} ) \Omega_{cdka} \Omega_{bacl} \Gamma_{ijbd} \right\} \\
%   IIb
\Gamma^{III_b}_{ijkl}  =& - \sum_{abcd} ( \bar{n}_b  {n}_c  {n}_d  +  {n}_b  \bar{n}_c  \bar{n}_d ) (1 - \hat P_{ij} ) (1 - \hat P_{kl} )  \nonumber\\
& \times \left( \Omega_{dcbk} \Omega_{biac} \Gamma_{jald} + 
 \Omega_{jcbd} \Omega_{balc} \Gamma_{diak} \right) \\
%    IIc
\Gamma^{III_c}_{ijkl}  =& - \frac{1}{2} \sum_{abcd} (   \bar{n}_a \bar{n}_b {n}_c + {n}_a  {n}_b  \bar{n}_c)  
 (1 - \hat P_{ij} )  (1 - \hat P_{kl} )  \nonumber\\
 & \times  \left(  \Omega_{abcl} \Omega_{idab} \Gamma_{cjkd} +  \Omega_{icab} \Omega_{abdl} \Gamma_{djkc}   \right)\\
%    IIIb
\Gamma^{IV_a}_{ijkl}  =&  -   \sum_{abcd} \left(   \bar{n}_c \bar{n}_d {n}_a  +  {n}_c  {n}_d  \bar{n}_a \right) \nonumber \\
& \times \left\{  (1\! -\! \hat P_{ij} ) \Omega_{aicd} \Omega_{dbkl} \Gamma_{jcba} \right. \nonumber \\
&+ \left.  (1 \!-\! \hat P_{kl} )  \Omega_{dcak} \Omega_{ijcb} \Gamma_{bald} \right\}\\
%    IIIa
\Gamma^{IV_b}_{ijkl}  =&  (1 - \hat P_{ij} ) (1 - \hat P_{kl} ) \sum_{abcd} (   \bar{n}_a {n}_b \bar{n}_c +  {n}_a  \bar{n}_b  {n}_c )   \nonumber\\
& \times \left(\Omega_{bica} \Omega_{jcld} \Gamma_{dabk}  + \Omega_{cabk} \Omega_{jdlc} \Gamma_{bida} \right)\\
%    IIIc
\Gamma^{IV_c}_{ijkl} =&   \frac{1}{2} (1 - \hat P_{ij} ) (1 - \hat P_{kl}) \sum_{abcd} (   \bar{n}_a \bar{n}_b {n}_d +{n}_a  {n}_b  \bar{n}_d ) \nonumber \\
&\times \left( \Omega_{abld} \Omega_{djck} \Gamma_{icab}  +\Omega_{idab} \Omega_{cjdk} \Gamma_{ablc}   \right)
\end{align}
$$

These may be expressed in a factorized form
$$
\begin{align}
\Gamma^{I}_{ijkl}  =& \sum_{a} \left\{(1\! -\! \hat P_{ij} )   \chi^{\epsilon}_{ia} \Gamma_{ajkl} % \right. \nonumber\\
+ (1\! -\! \hat P_{kl} )  \chi^{\epsilon}_{ak} \Gamma_{ijal}    \right\}\\
 %&+ \left. \left(1 - \hat P_{kl} \right)  \chi^{\epsilon}_{ak} \Gamma_{ijal}    \right\}\\
%
\Gamma^{II}_{ijkl}  =&\sum_{a}  \left\{
(1 \!-\! \hat P_{ij} )   \chi^{\zeta}_{aj}  \Omega_{iakl} % \right.\nonumber\\
%& - \left.
-(1 \!-\! \hat P_{kl} )  \chi^{\zeta}_{ak}  \Omega_{ijal}  \right\}\\
%
\Gamma^{III_a}_{ijkl}  =& 
-\!\sum_{ab} \bigl\{ (1\! -\! \hat P_{ij} ) \chi^{\eta}_{ijab}  \Gamma_{abkl} %\nonumber\\
  + (1 \!- \!\hat P_{kl} )  \Gamma_{ijab}  \chi^{\eta}_{abkl} \bigr\} \\
%-\left(1 - \hat P_{ij} \right) \sum_{ab}\chi^{\eta}_{ijab}  \Gamma_{abkl} \nonumber\\
%& - \left(1 - \hat P_{kl} \right) \sum_{ab}  \Gamma_{ijab}  \chi^{\eta}_{abkl}  \\
%
\Gamma^{III_b}_{ijkl} =& - (1\! -\! \hat P_{ij}) (1\! -\! \hat P_{kl} ) % \nonumber\\
%&\times
\sum_{ab} \left( \chi^{\eta}_{bkai} \Gamma_{jbla} + \chi^{\eta}_{lajb}  \Gamma_{aibk}  \right)\\
%
\Gamma^{III_c}_{ijkl}  =& - \frac{1}{2} \left(1 - \hat P_{ij} \right)  (1\! - \!\hat P_{kl} ) \sum_{ab}  \chi^{\theta}_{iabl} \Gamma_{bjka}   \\
%
\Gamma^{IV_a}_{ijkl}  =&  -\!\sum_{ab} \left\{  (1\!- \!\hat P_{ij} ) \chi^{\kappa}_{ijab} \Omega_{bakl}  \right. %\nonumber\\ 
%&+ 
-\left.  (1\! -\! \hat P_{kl} ) \Omega_{ijab} \chi^{\kappa}_{klba} \right\}  \\
%
\Gamma^{IV_b}_{ijkl}  =&  (1 \!- \!\hat P_{ij} )  (1 \!-\! \hat P_{kl} ) %\nonumber \\
%&\times 
\sum_{ab} \left( \chi^{\iota}_{aibk} \Omega_{jbla} - \chi^{\iota}_{akbi} \Omega_{jalb}  \right ) \\
%
%\Gamma^{IV_c}_{ijkl}   =&  - \frac{1}{2} \left(1 - \hat P_{ij} \right) \left(1 - \hat P_{kl} \right)  \sum_{ab} \chi^{\lambda}_{jalb}  \Omega_{biak} 
\Gamma^{IV_c}_{ijkl}   =&   \frac{1}{2} (1\! -\! \hat P_{ij} ) 
 (1 \!-\! \hat P_{kl} )  \sum_{ab} \chi^{\lambda}_{ialb}  \Omega_{bjak} 
\end{align}
$$
where we have defined additional intermediates
$$
\begin{align}
\chi^{\epsilon}_{ij} =& \frac{1}{2} \sum_{abc} \left( \bar{n}_a \bar{n}_b {n}_c    +  {n}_a  {n}_b \bar{n}_c \right)  \Omega_{ciab} \Omega_{abcj} \\  % I   ab
%
\chi^{\zeta}_{ij}  =&\frac{1}{2} \sum_{abc} \left(    {n}_a {n}_b \bar{n}_c + \bar{n}_a  \bar{n}_b {n}_c \right) \Gamma_{aibc}  \Omega_{bcaj} \\   % IV ab
%
\chi^{\eta}_{ijkl} =& \sum_{ab} \left( \bar{n}_a  {n}_b  \bar{n}_k  +  {n}_a  \bar{n}_b  {n}_k \right)  \Omega_{iabl} \Omega_{bjka} \\    % II  ac   % II  bd
%
\chi^{\theta}_{ijkl}  =&  \sum_{ab} \left(   {n}_a {n}_b \bar{n}_k  +  \bar{n}_a  \bar{n}_b {n}_k \right. \nonumber\\
& ~~~+ \left.   {n}_a {n}_b \bar{n}_j + \bar{n}_a  \bar{n}_b {n}_j  \right)   \Omega_{ijab} \Omega_{abkl}  \\   % II  ef
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chi^{\iota}_{ijkl} =& \sum_{ab} \left(   \bar{n}_a {n}_b \bar{n}_k  +  {n}_a  \bar{n}_b  {n}_k \right)  \Omega_{bika} \Gamma_{iabl} \\     % III ab
%
\chi^{\kappa}_{ijkl}  =&\sum_{ab} \left(   \bar{n}_a {n}_b {n}_k  +  {n}_a  \bar{n}_b  \bar{n}_k \right)  %\nonumber \\
\Omega_{ajbl} \Gamma_{ibka} \\     % III cd
%
\chi^{\lambda}_{ijkl}  =& \sum_{ab} \bigl\{ \left(   \bar{n}_a \bar{n}_b {n}_l +  {n}_a  {n}_b  \bar{n}_l \right)  \Gamma_{ijab} \Omega_{abkl}  \nonumber\\
                        &~~+   \left(   \bar{n}_a \bar{n}_b {n}_j +  {n}_a  {n}_b  \bar{n}_j \right)  \Omega_{ijab} \Gamma_{abkl} \bigr\}
% III ef
\end{align}
$$
We note that the intermediates $\chi$ do not in general have Hermitian symmetry or permutation symmetry under exchange of indices.


\section{J-coupled factorization version of double commutators}
To take advantage of rotational symmetry, we work with J-coupled matrix elements.
We work with un-normalized coupled states so that
$$
\begin{equation}
    A_{ijkl}^{JM} = \sum_{m_im_jm_im_l} \mathcal{C}^{JM}_{j_im_i,j_jm_j}
     \mathcal{C}^{JM}_{j_km_k,j_lm_l}
    A_{ijkl}
\end{equation}
$$
where the $\mathcal{C}$ are Clebsch-Gordan coefficients, and the uncoupled matrix element $A_{ijkl}$ depends on the projections $m_i,m_j,m_k,m_l$.
For the rotationally invariant operators considered in this work, $A^{JM}_{ijkl}$ is independent of the total projection $M$ and so $M$ is not explicitly indicated.
These expressions were obtained with the aid of the \texttt{amc} code~\cite{Tichai2020}.
In several cases, the expressions are simplified and made amenable to matrix multiplication if we use different coupling schemes.
In addition to the standard coupling scheme, we employ both Pandya-transformed and cross-coupled matrix elements, illustrated in Fig.~\ref{fig_coupling}.
Pandya-transformed matrix elements are indicated with a single bar
$$
\begin{align} \label{Pandya}
\bar A^J_{i \bar l k \bar j} = - \sum_{J^\prime} (\hat{J}^\prime)^2 \begin{Bmatrix} j_i & j_l & J \\ j_k & j_j & J^\prime \end{Bmatrix}  A^{J^\prime}_{i j k l} 
\end{align}
$$
where the braces indicate the six-J coefficient,
$j_a$ denotes the angular momentum of the orbit $a$
and we use the usual notation $\hat J \equiv \sqrt{2J+1}$.
Cross-coupled matrix elements are indicated with a double bar
$$
\begin{align} \label{cross-coupled_Pandya}
 \overline{\overline{A}}^J_{j \bar l k \bar i} &= \sum_{J^\prime} (\hat{J}^\prime)^2 \begin{Bmatrix} j_j & j_l & J \\ j_k & j_i & J^\prime \end{Bmatrix}  (-1)^{j_i +j_j -J^\prime}   A^{J^\prime}_{i j k l} \nonumber \\
 \overline{\overline{A}}^J_{i \bar k l \bar j} &= \sum_{J^\prime} (\hat{J}^\prime)^2 \begin{Bmatrix} j_i & j_k & J \\ j_l & j_j & J^\prime \end{Bmatrix}  (-1)^{j_k +j_l -J^\prime}   A^{J^\prime}_{i j k l} 
\end{align}
$$
In these definitions, we have not assumed Hermitian or permutation symmetries in the matrix elements $A$.


\begin{figure}[ht]
\centering
\input{Diagrams/coupling_diagrams.tikz}
\caption{\label{fig_coupling}Coupling schemes used in $J$-coupled expressions. (a)~Standard coupling, (b)~Pandya-transformed matrix elements, (c)~cross-coupled matrix elements.}
\end{figure}

The one-body matrix elements are
$$
\begin{align}  
 f^I_{ij}  =& \delta_{ j_i j_j} \hat{{j}}_{i}^{-2}  \sum_{ab J}  \hat{{J}}^{2}  \*  \chi^{\alpha}_{ab}  \Gamma^{J}_{biaj}      \\
    %
 f^{II}_{ij}  =& \delta_{ j_i j_j} ~\hat j_i^{-2}  \sum_{ab J}  \hat J^2 \left(\chi^{\beta}_{ab}  -   \chi^{\beta}_{ba}  \right) \Omega^{J}_{biaj}   \\
 %
 f^{III_a}_{ij}  =&  \delta_{ j_i j_j}~\hat j_i^{-2} \sum_{abc J} \left(  \bar \chi^{\gamma~J}_{i \bar c a \bar b} \Gamma^{J}_{a \bar b j \bar c} - \bar \chi^{\gamma~J}_{c \bar j a \bar b} \Gamma^{J}_{a \bar b c \bar i} \right)  \\
  %
 f^{III_b}_{ij}  =& \delta_{ j_i j_j} ~\hat j_i^{-2} \sum_{abc J} \left(\chi^{\delta~J}_{ciab}  \Gamma^{J}_{abcj} - \chi^{\delta~J}_{abcj} \Gamma^{J}_{ciab} \right)
\end{align}
$$
with intermediates
$$
\begin{align}  
 &\chi^{\alpha}_{ij} = \!  \sum_{abc J} \frac{\hat{J}^2}{2\hat{j}_{i}^{2}}   \left( \bar{n}_a  \bar{n}_b  {n}_c {n}_i  -  {n}_a  {n}_b  \bar{n}_c  \bar{n}_i \right) %\nonumber\\
  \Omega^{J}_{ciab} \Omega^{J}_{abcj} \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%&\chi^\beta_{ij} =\! \sum_{abc J} \frac{\hat{J}^2}{2\hat{j}_{i}^{2}}   \left( \bar{n}_a  {n}_b  \bar{n}_c {n}_i  -  {n}_a  \bar{n}_b  {n}_c  \bar{n}_i \right) 
&\chi^\beta_{ij} =\! \sum_{abc J} \frac{\hat{J}^2}{2\hat{j}_{i}^{2}}   \left( \bar{n}_a  \bar{n}_b  {n}_c {n}_i  -  {n}_a  {n}_b  \bar{n}_c  \bar{n}_i \right) 
 \Omega^{J}_{ciab} \Gamma^{J}_{abcj} \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
&\bar{\chi}^{\gamma ~ J}_{i \bar j k \bar l} = \!  \sum_{ab}\hat J^2 \left( {n}_a \bar{n}_b \bar{n}_k  {n}_l    -  \bar{n}_a  {n}_b {n}_k  \bar{n}_l  \right)   \bar{\Omega}^{J}_{i \bar j  a \bar b} \bar{\Omega}^{J}_{a \bar b  k \bar l} \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
&\chi^{\delta~J}_{ijkl} =\!  \sum_{ab}\frac{\hat J^2}{4}\! \left(   {n}_a  {n}_b \bar{n}_k \bar{n}_l -    \bar{n}_a \bar{n}_b {n}_k {n}_l   \right) \Omega^{ J}_{ijab} \Omega^{ J}_{abkl}
\end{align}
$$

The $J$-coupled factorized expressions for the two-body matrix elements of the double commutator are
$$
\begin{align}
\Gamma^{I ~ J}_{ijkl}  =&  \sum_{a} \left\{ (1\! -\! \hat P^{J}_{ij} )   \chi^{\epsilon}_{ai} \Gamma^{J}_{ajkl} % \nonumber \\
%&
+ (1\! -\! \hat P^{J}_{kl} )   \chi^{\epsilon}_{ak} \Gamma^{J}_{ijal} \right\} \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\Gamma^{II ~ J}_{ijkl}  =&  \sum_{a}
\left\{ (1 \!-\! \hat P^{J}_{ij} )    \chi^{\zeta}_{aj}  \Omega^{J}_{iakl} 
 %\nonumber \\
%&
-(1\! -\! \hat P^{J}_{kl} )   \chi^{\zeta}_{ak}  \Omega^{J}_{ijal} \right\}
\\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\Gamma^{III_a J}_{ijkl}  =&   %\nonumber\\
%&\times
\sum_{ab} \left\{  (1 - \hat P^{J}_{ij} )  \chi^{\eta~J}_{ijab}  \Gamma^{J}_{abkl} \right. %\nonumber \\ 
%&
+  \left. (1 - \hat P^{J}_{kl} ) \Gamma^{J}_{ijab} \chi^{\eta~J}_{a b k l} \right\}\\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\overline{\overline{ \Gamma}}^{III_bJ}_{j \bar l k \bar i} =& (1 - \hat P^{J}_{ij} )  (1 - \hat P^{J}_{kl} ) % \nonumber \\
%&\times
\sum_{ab} \overline{\overline{ \Gamma}}^{J}_{j \bar l a \bar b}  \left(   \overline{\overline{\chi}}^{\eta~J}_{ a \bar b k \bar i} 
  + \overline{\overline{ \chi}}^{\eta~J}_{k \bar i  b \bar a}  \right) \\ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bar \Gamma^{III_c J}_{i \bar l k \bar j}  =& 
\frac{1}{2} (1 - \hat P^{J}_{ij} )  (1 - \hat P^{J}_{kl} ) \sum_{ab}  \bar \chi^{\theta~J}_{i\bar l a \bar b}  \bar \Gamma^{J}_{a \bar b k \bar j}
\\
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\Gamma^{IV_a J}_{ijkl}  =&
-\sum_{ab} \left\{
(1 - \hat P^{J}_{ij} ) \chi^{\kappa~J}_{ijab} ~ \Omega^{J}_{abkl}    \right.  \nonumber \\
&\hspace{2em}+  \left. (1 - \hat P^{J}_{kl} ) \Omega^{J}_{ijab} ~ \chi^{\kappa~J}_{klab}  \right\} \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\overline{\overline{ \Gamma}}^{IV_b J}_{j \bar l k \bar i} =& (1 - \hat P^J_{ij} )  (1 - \hat P^J_{kl} ) % \nonumber \\ 
%& \times
\sum_{ab} \overline{\overline{  \Omega}}^{J}_{j \bar l a \bar b} \left(   ~ \bar \chi^{\iota~J}_{a \bar b k \bar i}  
  -  \bar \chi^{\iota~J}_{ k \bar i a \bar b } \right) \\ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bar \Gamma^{IV_c~J}_{i \bar l k \bar j}  =&   \frac{1}{2} (1 - \hat P^J_{ij} )  (1 - \hat P^J_{kl} ) \sum_{ab} \bar \chi^{\lambda~J}_{i \bar l a \bar b}  ~ \bar\Omega^{J}_{a\bar k j \bar  b}  
\end{align}
$$
and the $J$-coupled intermediates are
$$
\begin{align}
\chi^{\epsilon}_{ij} =& \frac{1}{2\hat{j}^2_j} \sum_{abc J} \hat{J}^2  \left( \bar{n}_a \bar{n}_b {n}_c    +  {n}_a  {n}_b \bar{n}_c \right)  \Omega^{J}_{ciab} \Omega^{J}_{abcj} \\  % I   ab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chi^{\zeta}_{ij}  =&\frac{1}{2} \sum_{abc J} \hat{J}^2 j^{-2}_j \left( \bar{n}_a  \bar{n}_b {n}_c  + {n}_a {n}_b \bar{n}_c  \right) \nonumber \\
 &\times  \Gamma^{J}_{aibc}  \Omega^{J}_{bcaj} \\   % IV ab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bar \chi^{\eta~J}_{i\bar jk \bar l} =& \sum_{ab} \left( \bar{n}_a  {n}_b  \bar{n}_k  +  {n}_a  \bar{n}_b  {n}_k \right) \bar \Omega^{J}_{i \bar j a \bar b} \bar \Omega^{J}_{a \bar b k \bar l} \\    % II  ac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chi^{\theta~J}_{ijkl}  =&  \sum_{ab} \left(   {n}_a {n}_b \bar{n}_k  +  \bar{n}_a  \bar{n}_b {n}_k + {n}_a {n}_b \bar{n}_j  \right. \nonumber\\
&+  \left. \bar{n}_a  \bar{n}_b {n}_j \right)  \Omega^{J}_{ijab} \Omega^{J}_{abkl}  \\ 
% II  e f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bar \chi^{\iota~J}_{i\bar j k \bar l} =& \sum_{ab} \left(   \bar{n}_a {n}_b \bar{n}_k  +  {n}_a  \bar{n}_b  {n}_k \right) \bar\Gamma^{J}_{i \bar j a \bar b}   \bar \Omega^{J}_{a \bar b k \bar l} \\     % III ab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\overline{\overline{  \chi}}^{\kappa~J}_{i \bar jk \bar l}  =&   \sum_{ab} \left(   \bar{n}_a {n}_b {n}_l  +  {n}_a  \bar{n}_b  \bar{n}_l \right) 
    \overline{\overline{  \Omega}}^{J}_{i \bar j b \bar a} \overline{\overline{  \Gamma}}^{J}_{b \bar a k \bar l} \\     % III cd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\chi^{\lambda~}_{ijkl}  =& \sum_{ab} \left(   \bar{n}_a \bar{n}_b {n}_l +  {n}_a  {n}_b  \bar{n}_l \right)  \Gamma^{J}_{ijab} \Omega^{J}_{abkl}    \nonumber\\   
 &+\sum_{ab} \left(   \bar{n}_a \bar{n}_b {n}_j +  {n}_a  {n}_b  \bar{n}_j \right)  \Omega^{J}_{ijab} \Gamma^{J}_{abkl}  
% III ef
\end{align}
$$
where $\hat P^{J}_{ij} \equiv (-1)^{ j_i  + j_j - J  } P^{J}_{ij} $ is the J-coupled permutation operator.


We note that the intermediate operators appearing in 
equations~\eqref{Jcoupled_Factorized_DoubleCommutator_twobody} do not always have the same coupling scheme with which they are defined in equations~\eqref{chi_2b_Jcoupled}.
The translation between these can be obtained from equations~\eqref{Pandya} and \eqref{cross-coupled_Pandya}.

In practice, one can derive an equation to manage the direct transformation instead of taking two steps.