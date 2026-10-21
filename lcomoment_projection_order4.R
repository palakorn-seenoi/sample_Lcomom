# ================================================================
# Projection (first-projection, Hajek) covariance estimator Omega_hat_n
# for directional sample L-comoments of orders up to four.
#
# Method: the first projection of each U-statistic is estimated at every
# observation by the leave-one-out identity
#
#   M_hat_{m,1}(Z_i) = { C(n,m) U_n - C(n-1,m) U_{n-1}^{(-i)} } / C(n-1,m-1),
#
# where U_{n-1}^{(-i)} is the same statistic computed without
# observation i. All leave-one-out statistics follow from one sort and
# prefix sums, so the cost is O(n log n) for the moments and O(n p^2)
# for the p x p matrix. The estimator Omega_hat_n is positive
# semidefinite by construction; Omega_hat_n / n estimates the leading
# term of the covariance of the moment vector.
#
# Moment vector (10): mean_X, mean_Y, L2_XY, L2_YX, L3_XY, L3_YX,
#                     L4_XY, L4_YX, L2_X, L2_Y
# (L_m_XY: X ranked by Y; L2_X, L2_Y: marginal L-scales.)
#
# Usage:
#   fit <- lcomoment_projection_cov4(x, y)
#   fit$moments; fit$Omega; fit$cov  (= Omega / n)
# No external packages are required.
# ================================================================

lweights <- function(n, m) {
  q <- seq_len(n)
  b <- list(rep(1 / n, n))
  for (r in 1:3) {
    t <- rep(1, n)
    for (j in 1:r) t <- t * (q - j) / (n - j)
    b[[r + 1]] <- t / n
  }
  switch(m,
         b[[1]],
         2 * b[[2]] - b[[1]],
         6 * b[[3]] - 6 * b[[2]] + b[[1]],
         20 * b[[4]] - 30 * b[[3]] + 12 * b[[2]] - b[[1]])
}

# v_sorted: values ordered by the ranking variable.
# Returns the full-sample statistic and the n leave-one-out statistics
# (in the same sorted order).
lstat_loo <- function(v_sorted, m) {
  n <- length(v_sorted)
  full <- sum(lweights(n, m) * v_sorted)
  w1 <- lweights(n - 1, m)
  a <- c(w1 * v_sorted[-n], 0)            # position p keeps index p   (p < q)
  b <- c(0, w1 * v_sorted[-1])            # position p moves to p - 1  (p > q)
  A <- c(0, cumsum(a))                    # A[q] = sum_{p < q} a_p
  B <- c(rev(cumsum(rev(b))), 0)          # B[q] = sum_{p >= q} b_p
  loo <- A[1:n] + B[2:(n + 1)]
  list(full = full, loo = loo)
}

lcomoment_projection_cov4 <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y); n <- length(x)
  if (n < 5) stop("n >= 5 is required for order-four projections")
  oy <- order(y); ox <- order(x)
  spec <- list(list(x, NULL, 1), list(y, NULL, 1),
               list(x, oy, 2), list(y, ox, 2),
               list(x, oy, 3), list(y, ox, 3),
               list(x, oy, 4), list(y, ox, 4),
               list(x, ox, 2), list(y, oy, 2))
  ord <- c(1, 1, 2, 2, 3, 3, 4, 4, 2, 2)
  U <- numeric(10); Mhat <- matrix(NA_real_, n, 10)
  for (k in seq_along(spec)) {
    v <- spec[[k]][[1]]; o <- spec[[k]][[2]]; m <- spec[[k]][[3]]
    if (m == 1) { U[k] <- mean(v); Mhat[, k] <- v; next }
    r <- lstat_loo(v[o], m)
    loo <- numeric(n); loo[o] <- r$loo
    U[k] <- r$full
    Mhat[, k] <- (choose(n, m) * r$full - choose(n - 1, m) * loo) / choose(n - 1, m - 1)
  }
  phi <- sweep(sweep(Mhat, 2, U), 2, ord, `*`)
  Omega <- crossprod(phi) / n
  nm <- c("mean_X", "mean_Y", "L2_XY", "L2_YX", "L3_XY", "L3_YX",
          "L4_XY", "L4_YX", "L2_X", "L2_Y")
  names(U) <- nm; dimnames(Omega) <- list(nm, nm)
  list(moments = U, Omega = Omega, cov = Omega / n)
}
