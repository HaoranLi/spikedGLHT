\name{spiked_glht}
\alias{spiked_glht}
\alias{print.spiked_glht}
\title{General Linear Hypothesis Tests with Spiked Noise Covariance}
\description{Tests \eqn{H_0: BC=D} in the model \eqn{Y=BX+E} using the
regularized LR, LH, and BNP statistics of Li et al. (2024).}
\usage{
spiked_glht(Y, X, C, D = NULL, rho = 0.5, K = NULL,
  calibration = c("normal", "bootstrap"), alpha = 0.05,
  B = 1999L, method = c("KN", "PLY"), conf_level = 0.95,
  delta = 0.05, K_max = NULL, adjustment = c("paper", "legacy"))
\method{print}{spiked_glht}(x, ...)
}
\arguments{
\item{Y}{Finite numeric response matrix, p rows by N subjects.}
\item{X}{Finite numeric design matrix, m rows by N subjects, of full row rank.
Include an intercept explicitly if required.}
\item{C}{Finite numeric contrast matrix, m rows by q columns, of full column rank.}
\item{D}{Finite numeric p by q null matrix. NULL denotes zero.}
\item{rho}{Maximum fraction of signal in the nonspiked subspace, in (0,1].
Corresponds to the paper's varrho and the scripts' c_0.}
\item{K}{Number of spikes. NULL selects it from residual eigenvalues.
Must be an integer between zero and min(p-1,N-m-1).}
\item{calibration}{Normal asymptotic approximation or parametric bootstrap.}
\item{alpha}{Significance level, strictly between zero and one.}
\item{B}{Positive integer number of bootstrap replicates. Used only for bootstrap.}
\item{method}{Spike selection: Kritchman--Nadler (KN) or Passemier--Li--Yao (PLY).
KN requires the optional RMTstat package.}
\item{conf_level}{Tracy--Widom quantile probability for KN, in (0,1).}
\item{delta}{Positive PLY penalty exponent.}
\item{K_max}{Maximum candidate spike count. Defaults to min(p-1,N-m-1,
max(1,floor(0.1*p))). Used for automatic spike selection.}
\item{adjustment}{Paper Algorithm 1, or legacy scripts' N/n eigenvalue inflation
and detection-boundary offset of 0.0001.}
\item{x}{A spiked_glht result.}
\item{...}{Arguments passed to print.data.frame.}
}
\details{
The model assumes independent Gaussian errors with a spiked covariance,
fixed predictor and contrast dimensions and a fixed number of separated,
supercritical spikes as p and N grow. Finite-sample calibration is approximate.
The bootstrap implements Theorem 3.2 and requires p >= 2*K.
It holds the estimated covariance and weights fixed across replicates.
Set the random seed with set.seed for reproducibility; the function uses R's
current random number stream.

QR residualization and an SVD avoid constructing the p by p covariance or
N by N projection matrix. Minimax weights solve Algorithm 2 using quadprog.
Zero eigenvector overlaps at the detection boundary are combined exactly
in the quadratic program. When K=0, the method uses the isotropic extension
with weight 1 and prior 1; the rho restriction is then inapplicable.

KN follows the supplied scripts' sequential threshold and upper search bound,
but allows zero spikes rather than forcing a minimum of one. PLY forwards both
delta and K_max. See the package's METHODS.md for differences from the scripts.

All three tests use the upper tail. Normal p-values use the standardized
statistics. Bootstrap p-values use (1 + number of replicates at least as
large as the observed raw statistic)/(B+1). Rejection uses the empirical
(1-alpha) raw-statistic quantile (R quantile type 7), as in the paper.
Consequently bootstrap rejection and p.value < alpha can differ near a
Monte Carlo boundary. No adjustment for reporting three tests is applied.
}
\value{An object of class spiked_glht with named vectors statistic (standardized),
p.value (chosen calibration), reject, raw.statistic, critical.value (raw scale),
and normal.p.value; weights and least.favorable.prior ordered as component0,
then the spikes; K, spikes, sigma_sq, theta1, theta2, alpha, calibration,
adjustment, rho, dimensions and call. bootstrap.statistics contains a B by 3
matrix of raw bootstrap statistics, or NULL for normal calibration.}
\references{Li, H., Aue, A., Paul, D., and Peng, J. (2024).
Testing General Linear Hypotheses Under a High-Dimensional Multivariate
Regression Model with Spiked Noise Covariance.
Journal of the American Statistical Association, 119(548), 2799--2810.
\doi{10.1080/01621459.2023.2278825}.}
\examples{
set.seed(7)
p <- 30; N <- 60
X <- rbind(rep(c(1, 0), each = N/2), rep(c(0, 1), each = N/2))
C <- matrix(c(1, -1), 2, 1)
Y <- matrix(rnorm(p * N), p, N)
Y[1, ] <- 4 * Y[1, ]
spiked_glht(Y, X, C, K = 1)
}
