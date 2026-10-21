# ================================================================
# Simulation experiments for the covariance estimators of sample
# L-comoments.
#
#   source("additional_experiments.R", chdir = TRUE)
#
# Requires lcomoment_fast_exact_cov.R. No external packages.
#
# Methods compared
#   * exactly unbiased covariance estimator (rank-prefix algorithm),
#   * its PSD projection,
#   * projection (first-projection) estimator Omega_hat_n / n,
# against a high-precision Monte Carlo reference covariance computed
# from R_REF replicated moment vectors. Under independence, the
# L-comoment block is also compared with its exact closed form.
# Also reported: first-order (delta-method) variance of L-comoment
# ratios relative to their Monte Carlo variance, coverage of Wald
# intervals for the ratios, and the elapsed time of the fast algorithm.
#
# Designs: IND    X ~ Exp(1), Y ~ Gumbel(0,1), independent
#          GAUSS  bivariate normal, rho = 0.6
#          GEVCOP Gaussian copula rho = 0.6, GEV(0,1,xi=0.2) margins
#
# Run time is dominated by the B covariance evaluations (about 0.26 s
# each at n = 50 on a laptop). B = 20000 takes many hours; B = 2000
# gives the same qualitative picture. The R_REF replications of the
# moment vector are cheap.
# ================================================================
source("lcomoment_fast_exact_cov.R")

DESIGNS <- c("IND", "GAUSS", "GEVCOP")
N_VALUES <- c(20L, 30L, 50L)
R_REF <- 1e6
B <- 2000L            # e.g. 20000 for higher precision
SEED <- 20261018L
RHO <- 0.6
XI <- 0.2
NUM <- 3:6; DEN <- c(7L, 8L, 7L, 8L)
ORD <- c(1, 1, 2, 2, 3, 3, 2, 2)
nm <- c("mean_X", "mean_Y", "L2_X_byY", "L2_Y_byX",
        "L3_X_byY", "L3_Y_byX", "L2_X", "L2_Y")

gev_from_z <- function(z, xi = XI) {
  # GEV(0,1,xi) quantile at pnorm(z), computed via log pnorm for stability
  ((-pnorm(z, log.p = TRUE))^(-xi) - 1) / xi
}

sample_design <- function(design, R, n) {
  if (design == "IND")
    return(list(x = matrix(rexp(R * n), R, n),
                y = matrix(-log(rexp(R * n)), R, n)))   # Gumbel(0,1)
  z1 <- matrix(rnorm(R * n), R, n); z2 <- matrix(rnorm(R * n), R, n)
  zy <- z1; zx <- RHO * z1 + sqrt(1 - RHO^2) * z2
  if (design == "GAUSS") return(list(x = zx, y = zy))
  list(x = gev_from_z(zx), y = gev_from_z(zy))
}

population_ratios <- function(design) {
  if (design == "IND") return(rep(0, 4))
  if (design == "GAUSS") return(c(RHO, RHO, 0, 0))
  # GEVCOP: Monte Carlo integration with 4e6 draws (numerical quadrature
  # gives tau2 = 0.6175, tau3 = 0.1179)
  zy <- rnorm(4e6); zx <- RHO * zy + sqrt(1 - RHO^2) * rnorm(4e6)
  x <- gev_from_z(zx); uy <- pnorm(zy); ux <- pnorm(zx)
  l2x <- mean(x * (2 * ux - 1))
  t2 <- mean(x * (2 * uy - 1)) / l2x
  t3 <- mean(x * (6 * uy^2 - 6 * uy + 1)) / l2x
  c(t2, t2, t3, t3)
}

exact_indep_block <- function(n) {
  vx <- 1; vy <- pi^2 / 6
  l2x <- 1/2; l3x <- 1/6; l2y <- log(2); l3y <- log(9/8)
  an <- (n + 1) / (3 * n * (n - 1)); bn <- (n + 1) * (n + 2) / (5 * n * (n - 1) * (n - 2))
  d <- n - 1
  matrix(c(an * vx, l2x * l2y / d, 0, l3x * l2y / d,
           l2x * l2y / d, an * vy, l2x * l3y / d, 0,
           0, l2x * l3y / d, bn * vx, l3x * l3y / d,
           l3x * l2y / d, 0, l3x * l3y / d, bn * vy), 4, 4, byrow = TRUE)
}

# Sample moment vectors for many samples (rows), via concomitant weights
moments_batch <- function(X, Y) {
  n <- ncol(X); q <- seq_len(n)
  w2 <- (2 * q - n - 1) / (n * (n - 1))
  w3 <- (6 * q^2 - 6 * (n + 1) * q + (n + 1) * (n + 2)) / (n * (n - 1) * (n - 2))
  oy <- t(apply(Y, 1, order)); ox <- t(apply(X, 1, order))
  ri <- row(X)
  Xc <- matrix(X[cbind(c(ri), c(oy))], nrow(X))   # X ordered by Y
  Yc <- matrix(Y[cbind(c(ri), c(ox))], nrow(Y))
  Xs <- t(apply(X, 1, sort)); Ys <- t(apply(Y, 1, sort))
  cbind(rowMeans(X), rowMeans(Y), Xc %*% w2, Yc %*% w2,
        Xc %*% w3, Yc %*% w3, Xs %*% w2, Ys %*% w2)
}

reference <- function(design, n) {
  chunk <- 1e5; M <- NULL
  for (k in seq_len(R_REF / chunk)) {
    s <- sample_design(design, chunk, n)
    M <- rbind(M, moments_batch(s$x, s$y))
  }
  Cref <- cov(M); Mc <- sweep(M, 2, colMeans(M))
  se_ref <- sqrt(pmax(crossprod(Mc^2) / nrow(M) - Cref^2, 0) / nrow(M))
  tau <- M[, NUM] / M[, DEN]
  list(m = colMeans(M), C = Cref, se = se_ref, var_tau = apply(tau, 2, var))
}

jac <- function(m) {
  J <- matrix(0, 4, 8)
  for (q in 1:4) { J[q, NUM[q]] <- 1 / m[DEN[q]]; J[q, DEN[q]] <- -m[NUM[q]] / m[DEN[q]]^2 }
  J
}
psd_proj <- function(C) { e <- eigen(C, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, 0), nrow(C)) %*% t(e$vectors) }
projection_omega <- function(U, S1, n) {
  phi <- sweep(sweep(S1, 2, choose(n - 1, ORD - 1), "/"), 2, U) %*% diag(ORD)
  crossprod(phi) / n
}
std_fro <- function(A, Cref) { d <- sqrt(diag(Cref)); sqrt(sum((A / outer(d, d))^2)) }

run_one <- function(design, n) {
  ref <- reference(design, n); Cref <- ref$C
  tau0 <- population_ratios(design)
  J0 <- jac(ref$m); taylor_over_mc <- diag(J0 %*% Cref %*% t(J0)) / ref$var_tau
  Ch <- Cp <- Om <- array(NA_real_, c(8, 8, B)); minev <- numeric(B)
  cov_psd <- cov_proj <- numeric(4)
  for (b in seq_len(B)) {
    s <- sample_design(design, 1, n)
    fit <- lcomoment_exact_cov_fast(s$x[1, ], s$y[1, ], cov_only = TRUE)
    C <- unname(fit$unbiased_cov); U <- unname(fit$moments)
    minev[b] <- min(eigen(C, symmetric = TRUE, only.values = TRUE)$values)
    Ch[, , b] <- C; Cp[, , b] <- psd_proj(C)
    Om[, , b] <- projection_omega(U, fit$S1, n) / n
    J <- jac(U); t <- U[NUM] / U[DEN]
    se1 <- sqrt(pmax(diag(J %*% Cp[, , b] %*% t(J)), 0))
    se2 <- sqrt(pmax(diag(J %*% Om[, , b] %*% t(J)), 0))
    cov_psd <- cov_psd + (abs(t - tau0) <= 1.959964 * se1)
    cov_proj <- cov_proj + (abs(t - tau0) <= 1.959964 * se2)
  }
  iu <- upper.tri(Cref, diag = TRUE)
  summ <- function(A) {
    mn <- apply(A, c(1, 2), mean); sdv <- apply(A, c(1, 2), sd)
    z <- (mn - Cref)[iu] / sqrt(sdv[iu]^2 / B + ref$se[iu]^2)
    c(max_abs_z = max(abs(z)),
      rel_bias = std_fro(mn - Cref, Cref) / std_fro(Cref, Cref),
      rel_rmse = sqrt(mean(apply(A, 3, function(a) std_fro(a - Cref, Cref)^2))) /
        std_fro(Cref, Cref))
  }
  out <- data.frame(design = design, n = n, nonpsd = mean(minev < -1e-12),
                    t(summ(Ch)), t(summ(Cp)), t(summ(Om)))
  names(out)[4:12] <- paste0(rep(c("unb_", "psd_", "proj_"), each = 3),
                             c("max_abs_z", "rel_bias", "rel_rmse"))
  rat <- data.frame(design = design, n = n,
                    t(setNames(taylor_over_mc, paste0("taylor_over_mc_", 1:4))),
                    t(setNames(cov_psd / B, paste0("cover_psd_", 1:4))),
                    t(setNames(cov_proj / B, paste0("cover_proj_", 1:4))))
  if (design == "IND") {
    Ce <- exact_indep_block(n); bl <- 3:6
    mn <- apply(Ch[bl, bl, ], c(1, 2), mean); sdv <- apply(Ch[bl, bl, ], c(1, 2), sd)
    iu4 <- upper.tri(Ce, diag = TRUE)
    out$max_abs_z_vs_exact_indep <- max(abs(((mn - Ce) / (sdv / sqrt(B)))[iu4]))
    out$max_rel_ref_vs_exact_indep <- max(abs(((Cref[bl, bl] - Ce) /
                                           sqrt(outer(diag(Ce), diag(Ce))))[iu4]))
  }
  list(est = out, ratio = rat)
}

set.seed(SEED)
res <- list()
for (d in DESIGNS) for (n in N_VALUES) {
  t0 <- proc.time()[["elapsed"]]
  r <- run_one(d, n); res[[length(res) + 1]] <- r
  cat(sprintf("%s n=%d done (%.0f s)\n", d, n, proc.time()[["elapsed"]] - t0))
  print(r$est, row.names = FALSE, digits = 4)
}
est <- do.call(rbind, lapply(res, function(r) { r$est[setdiff(names(r$est),
          c("max_abs_z_vs_exact_indep", "max_rel_ref_vs_exact_indep"))] }))
write.csv(est, "additional_estimator_comparison.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(res, `[[`, "ratio")), "additional_ratio_results.csv", row.names = FALSE)
ind <- do.call(rbind, lapply(res, function(r) if (!is.null(r$est$max_abs_z_vs_exact_indep))
  r$est[c("design", "n", "max_abs_z_vs_exact_indep", "max_rel_ref_vs_exact_indep")]))
write.csv(ind, "additional_independent_exact_check.csv", row.names = FALSE)

# ---- Scaling of the fast algorithm (pure R loop; slower than NumPy) ----
TIMING_N <- c(100L, 200L, 500L, 1000L)
tim <- data.frame(n = TIMING_N, seconds = NA_real_)
for (k in seq_along(TIMING_N)) {
  n <- TIMING_N[k]; x <- rnorm(n); y <- RHO * x + sqrt(1 - RHO^2) * rnorm(n)
  tim$seconds[k] <- median(replicate(3, system.time(
    lcomoment_exact_cov_fast(x, y, cov_only = TRUE))[["elapsed"]]))
  print(tim[k, ])
}
tim$slope <- coef(lm(log(seconds) ~ log(n), data = tim))[2]
write.csv(tim, "additional_timing_R.csv", row.names = FALSE)
cat("Done. Files: additional_*.csv\n")
