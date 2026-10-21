# ================================================================
# Self-test for the covariance estimators of sample L-comoments.
#
#   source("self_test_fast_cov.R", chdir = TRUE)
#
# Test 1  Exactly unbiased covariance estimator:
#         the rank-prefix algorithm (lcomoment_exact_cov_fast) is
#         compared with direct enumeration of the disjoint-subset
#         estimator for n = 6, 7, 8.
# Test 2  Projection (first-projection) estimator:
#         Omega_hat_n from the leave-one-out algorithm
#         (lcomoment_projection_cov4) is compared with Omega_hat_n
#         formed from the incidence sums S1 returned by
#         lcomoment_exact_cov_fast, for the common moments.
# Test 3  Order-four leave-one-out projections are compared with
#         brute-force averages of the kernels containing each
#         observation.
#
# Direct enumeration costs O(n^6) and is used only for small n.
# No external packages are required.
# ================================================================
source("lcomoment_fast_exact_cov.R")
source("lcomoment_projection_order4.R")

# ---- Direct enumeration of the disjoint-subset estimator -------------
# For U-statistics U_m and U_s with kernels of orders m and s, the
# exactly unbiased estimator of Cov(U_m, U_s) is
#   U_m U_s - {C(n,m) C(n-m,s)}^{-1} sum_{I cap J = empty} h_m(I) h_s(J).
direct_cov <- function(x, y) {
  n <- length(x)
  H1 <- lapply(seq_len(n), function(i) list(idx = i, h = c(x[i], y[i])))
  pr <- combn(n, 2)
  H2 <- lapply(seq_len(ncol(pr)), function(k) {
    i <- pr[1, k]; j <- pr[2, k]
    dx <- x[j] - x[i]; dy <- y[j] - y[i]
    list(idx = c(i, j),
         h = c(dx * sign(dy) / 2, dy * sign(dx) / 2, abs(dx) / 2, abs(dy) / 2))
  })
  tr <- combn(n, 3)
  H3 <- lapply(seq_len(ncol(tr)), function(k) {
    t3 <- tr[, k]
    zx <- x[t3][order(y[t3])]          # X concomitants ordered by Y
    zy <- y[t3][order(x[t3])]          # Y concomitants ordered by X
    list(idx = t3,
         h = c((zx[1] - 2 * zx[2] + zx[3]) / 3, (zy[1] - 2 * zy[2] + zy[3]) / 3))
  })
  H <- list(H1, H2, H3)
  # positions in the moment vector of the kernels of each order
  slots <- list(c(1L, 2L), c(3L, 4L, 7L, 8L), c(5L, 6L))
  U <- numeric(8)
  for (m in 1:3)
    U[slots[[m]]] <- rowMeans(matrix(sapply(H[[m]], function(a) a$h),
                                     nrow = length(slots[[m]])))
  C <- matrix(0, 8, 8)
  for (m in 1:3) for (s in 1:3) {
    acc <- matrix(0, length(slots[[m]]), length(slots[[s]]))
    for (a in H[[m]]) for (b in H[[s]])
      if (!any(a$idx %in% b$idx)) acc <- acc + outer(a$h, b$h)
    C[slots[[m]], slots[[s]]] <-
      outer(U[slots[[m]]], U[slots[[s]]]) - acc / (choose(n, m) * choose(n - m, s))
  }
  list(moments = U, cov = (C + t(C)) / 2)
}

# ---- Projection estimator from the incidence sums S1 ------------------
projection_from_S1 <- function(fit) {
  n <- fit$n
  ord <- c(1, 1, 2, 2, 3, 3, 2, 2)
  U <- unname(fit$moments)
  phi <- sweep(sweep(fit$S1, 2, choose(n - 1, ord - 1), "/"), 2, U)
  phi <- sweep(phi, 2, ord, "*")
  crossprod(phi) / n
}

# ---- Order-four kernel average containing observation i ---------------
brute_projection4 <- function(x, y, i) {
  n <- length(x)
  others <- combn(setdiff(seq_len(n), i), 3)
  vals <- apply(others, 2, function(J) {
    I4 <- c(i, J)
    z <- x[I4][order(y[I4])]
    (-z[1] + 3 * z[2] - 3 * z[3] + z[4]) / 4
  })
  mean(vals)
}

tol <- 1e-10
ok <- TRUE
set.seed(1)

cat("Test 1: rank-prefix algorithm vs direct enumeration\n")
for (n in 6:8) {
  x <- rnorm(n); y <- 0.6 * x + rnorm(n)
  f <- lcomoment_exact_cov_fast(x, y, cov_only = TRUE)
  d <- direct_cov(x, y)
  e_mom <- max(abs(unname(f$moments) - d$moments))
  e_cov <- max(abs(unname(f$unbiased_cov) - d$cov))
  cat(sprintf("  n = %d: max |moment diff| = %.2e, max |cov diff| = %.2e\n",
              n, e_mom, e_cov))
  ok <- ok && e_mom < tol && e_cov < tol
}

cat("Test 2: leave-one-out projection vs incidence-sum projection\n")
common <- c(1:6, 9:10)      # positions in the order-four moment vector
for (n in c(8, 15, 40)) {
  x <- rexp(n); y <- x + rnorm(n)
  f <- lcomoment_exact_cov_fast(x, y, cov_only = TRUE)
  p4 <- lcomoment_projection_cov4(x, y)
  e_mom <- max(abs(unname(p4$moments[common]) - unname(f$moments)))
  e_om <- max(abs(unname(p4$Omega[common, common]) - projection_from_S1(f)))
  cat(sprintf("  n = %d: max |moment diff| = %.2e, max |Omega diff| = %.2e\n",
              n, e_mom, e_om))
  ok <- ok && e_mom < tol && e_om < tol
}

cat("Test 3: order-four projections vs brute-force kernel averages\n")
for (n in c(6, 9)) {
  x <- rnorm(n); y <- x + rnorm(n)
  oy <- order(y)
  r <- lstat_loo(x[oy], 4)
  loo <- numeric(n); loo[oy] <- r$loo
  Mhat <- (choose(n, 4) * r$full - choose(n - 1, 4) * loo) / choose(n - 1, 3)
  brute <- sapply(seq_len(n), function(i) brute_projection4(x, y, i))
  e4 <- max(abs(Mhat - brute))
  cat(sprintf("  n = %d: max |projection diff| = %.2e\n", n, e4))
  ok <- ok && e4 < tol
}

cat(if (ok) "All tests passed.\n" else "SOME TESTS FAILED.\n")
