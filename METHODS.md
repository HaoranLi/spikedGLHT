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

## Differences from the supplied scripts

- Algorithm 1 in the published paper solves the spike equation using tau/sigma².
  `equation_solver.R` additionally inflates it by N/n. The package makes this
  explicit with `adjustment = "paper"` (default) and `"legacy"`.
- Legacy mode preserves the 0.0001 lower offset used in the scripts. Paper mode
  uses the exact detection boundary in Algorithm 1's final step.
- At exactly zero overlap the prior QP has redundant columns. The package
  combines those columns, allocating their mass to an undetectable spike with
  pi_0=0. This preserves the objective and feasibility without a ridge penalty.
- KN retains the scripts' candidate-k noise estimate and p-k threshold, but
  returns zero when the first candidate fails. The scripts force at least one.
  The candidate search bound is always respected. This is a practical sequential
  implementation; supplementary threshold indexing is ambiguous at k=0.
- K=0 is an isotropic extension with lambda_0=1 and pi_0=1. Its prior is not in
  the rho<1 restricted family (which is empty at K=0); rho is explicitly ignored.
- PLY receives the actual delta and search bound; the original caller omitted them.
- The nonzero-null calculation uses A = YQ - D(C'(XX')^-1 C)^-1/2 and M=A'Omega A/p.
  The scripts subtract p by q D_C from q by p Q'Y', which is not conformable.
- Normal p-values are returned as upper-tail probabilities. The research
  `return_stat=FALSE` path returns rejection flags despite its p-value comment.
- Bootstrap calibration is optional and the replicate count is configurable;
  the real-data script always uses 5000 draws.
- LR uses sum(log1p(eigenvalues)) to avoid the scripts' product overflow.
- Projection matrices and dense p by p scaling matrices are avoided. SVD retains
  the full spectrum including zero eigenvalues when computing the bulk mean.
- Domain errors (eigenvalues of M at or below -1), rank failures, incompatible
  matrix dimensions and invalid tuning parameters produce explicit errors.

## Validation

`tests/legacy-fixture.rds` was generated from the unmodified simulation
`Proposed.R`, `equation_solver.R`, and `Minimax.R` using seed 20261018,
p=40, N=80, m=3, q=2, K=2, and c_0=0.5. It stores the inputs and original
statistics, weights, prior and covariance estimates. `tests/core.R` compares
all of these quantities with legacy mode to tolerance 1e-8.

Other checks independently evaluate paper Algorithm 1 and moment formulas,
rotation invariance, nonzero-null translation, bootstrap reproducibility and
its LH mean, zero and one spike, p>N, automatic spike selection, detection
boundaries, and input validation. These checks establish numerical consistency,
not a new proof of the asymptotic theory or reproduction of every paper table.
