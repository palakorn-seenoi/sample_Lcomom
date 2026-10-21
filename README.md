Exact finite-sample covariance of sample L-comoments: R code
Authors: Palakorn Seenoi and Jeong-Soo Park
This repository contains R code for estimating the finite-sample covariance matrix of sample L-comoments (directional L-comoments) and for the numerical experiments that assess these estimators.
Requirements: R (4.0 or later). No external packages are needed.
Put all files of the R folder in one folder and run each script with:
```r
source("<script name>", chdir = TRUE)
```
Moment vector
The covariance functions work with the vector
```text
mean(X), mean(Y), L2[X|rank Y], L2[Y|rank X],
L3[X|rank Y], L3[Y|rank X], L2(X), L2(Y),
```
that is, the sample means, the directional sample L-comoments of orders two and three, and the marginal sample L-scales, computed from i.i.d. pairs (X, Y) with continuous, tie-free margins. X and Y may be dependent. The projection estimator also covers order four.
Functions
`lcomoment_fast_exact_cov.R`
```r
lcomoment_exact_cov_fast(x, y, pd_floor_scaled = NULL, cov_only = FALSE)
```
Exactly unbiased covariance estimator (disjoint-subset U-statistic estimator), computed by the rank-prefix (inclusion--exclusion) algorithm in O(n^2) time.
PSD projection and positive-definite eigenvalue correction.
First-order (delta-method) covariance of the L-comoment ratios tau2 and tau3.
Incidence sums S1, from which the projection estimator is formed.
`lcomoment_projection_order4.R`
```r
lcomoment_projection_cov4(x, y)
```
Projection (first-projection, Hajek) covariance estimator Omega_hat_n for directional L-comoments of orders up to four, computed by a leave-one-out identity with prefix sums.
Example
```r
source("lcomoment_fast_exact_cov.R")
set.seed(17)
x <- rnorm(60); y <- 0.6 * x + rnorm(60)
fit <- lcomoment_exact_cov_fast(x, y)
fit$unbiased_cov          # exactly unbiased 8 x 8 covariance
fit$pd_corrected_cov      # positive-definite version
fit$ratio_estimates       # tau2 and tau3 ratios
fit$ratio_cov_delta_psd   # first-order covariance of the ratios
```
Correctness check
`self_test_fast_cov.R`
Compares (1) the rank-prefix algorithm with direct enumeration of the disjoint-subset estimator for n = 6, 7, 8, (2) the leave-one-out projection estimator with the projection estimator formed from S1, and (3) order-four leave-one-out projections with brute-force kernel averages. Prints "All tests passed." on success. Runs in seconds.
Numerical experiments
`benchmark_n20.R`, `benchmark_n30.R`, `benchmark_n50.R`
Computation time of the full 8 x 8 exactly unbiased covariance matrix: direct enumeration of the disjoint-subset estimator versus the rank-prefix algorithm, on the same bivariate Gaussian sample (rho = 0.6, seed 20261017). Also reports the largest absolute difference between the two matrices.
Output: `lcomoment_n<N>_*.csv`
`monte_carlo_lcomoment_cov.R`
Monte Carlo check of the exactly unbiased estimator for n = 20, 30, 50 (bivariate Gaussian, rho = 0.6, 1000 samples, seed 20261018): average of the covariance estimates versus the empirical covariance of the moment vectors, frequency of non-PSD estimates, and minimum eigenvalues.
Output: `n<N>_*.csv`, `mc_summary_all.csv`, `mc_entrywise_all.csv`
`additional_experiments.R`
Comparison of three covariance estimators (exactly unbiased, its PSD projection, and the projection estimator Omega_hat_n / n) with a high-precision Monte Carlo reference covariance, for three designs:
independent Exp(1) and Gumbel(0,1) margins,
bivariate Gaussian (rho = 0.6), and
a Gaussian copula with GEV(0,1,0.2) margins.
Under independence, the L-comoment block is also compared with its exact closed form. Also reports the accuracy of the first-order ratio variance, coverage of Wald intervals for the ratios, and the elapsed time of the rank-prefix algorithm for n up to 1000.
The default B = 2000 takes an hour or more on a laptop; B = 20000 gives higher precision but takes many hours.
Output: `additional_*.csv`
