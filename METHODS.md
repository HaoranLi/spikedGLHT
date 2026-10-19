# Method and implementation notes

Reference: Li, Aue, Paul and Peng (2024), JASA 119(548), 2799–2810,
DOI 10.1080/01621459.2023.2278825. The published PDF in the parent project
was used as the authoritative source; the LaTeX and research R scripts were
cross-checked for implementation details.

| Paper | Implementation |
| --- | --- |
| Section 2.2, equation (2.6) | QR residualization, contrast normalization, and low-rank evaluation of M |
| Algorithm 1, equation (3.12) | `.adjusted()` in `R/covariance.R` |
| Theorem 3.1, eigenvector overlap | `.weight_geometry()` in `R/weights.R` |
| Corollary 3.1, moments | `theta1`, `theta2` in `R/test.R` |
| Corollary 3.2, equation (3.7) | Normalized LR, LH and BNP in `R/test.R` |
| Theorem 3.2, equations (3.5)–(3.6) | `.bootstrap()` in `R/test.R` |
| Proposition 4.1 and Algorithm 2 | `.minimax()` in `R/weights.R` |
| Supplementary Section S.3 | KN and PLY selection in `.select_K()` |

## Implementation details

- Spike estimation follows Algorithm 1 using tau/sigma² and the exact detection
  boundary in its final step.
- At exactly zero overlap the prior QP has redundant columns. The package
  combines those columns, allocating their mass to an undetectable spike with
  pi_0=0. This preserves the objective and feasibility without a ridge penalty.
- KN uses a bounded sequential Tracy–Widom threshold search and permits zero
  spikes. PLY uses the specified delta and maximum candidate count.
- K=0 is an isotropic extension with lambda_0=1 and pi_0=1. Its prior is not in
  the rho<1 restricted family (which is empty at K=0); rho is ignored.
- For a nonzero null, A = YQ - D(C'(XX')^-1 C)^-1/2 and M=A'Omega A/p.
- Normal p-values are upper-tail probabilities. Bootstrap calibration implements
  Theorem 3.2 with a configurable replicate count.
- LR uses sum(log1p(eigenvalues)) for numerical stability.
- QR residualization and SVD avoid dense projection and scaling matrices. The
  bulk mean includes the zero residual eigenvalues when p exceeds sample size.

## Validation

The tests generate deterministic synthetic data and independently evaluate
Algorithm 1 and Corollary 3.1. They also check rotation invariance, nonzero-null
translation, bootstrap reproducibility and LH moments, zero and one spike,
p>N, automatic spike selection, detection boundaries, and input validation.
These checks establish numerical consistency, not a reproduction of every
paper table.
