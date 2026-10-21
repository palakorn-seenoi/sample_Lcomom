# Fast exactly unbiased finite-sample covariance for sample L-comoments.
# 17 October 2026. Research implementation for JMVA L-comoment manuscript.
# 18 October 2026: also returns S1 (n x 8 incidence sums), used for the
#   projection estimator Omega_hat_n of Theorem 4.3. No other change.
#
# Target vector (in order):
# [mean(X), mean(Y), L2[X|rank Y], L2[Y|rank X],
#  L3[X|rank Y], L3[Y|rank X], L2(X), L2(Y)].
# All 8x8 entries are exactly unbiased estimators of finite-sample covariance.
# Algorithm uses inclusion--exclusion of index overlaps, NOT independent-pair enumeration.
# Runtime O(n^2); memory O(n^2) from bivariate rank-prefix sums.
# Assumptions: iid (X,Y) blocks, continuous tie-free margins, n>=6.
# Margins may be dependent within a block.
# No Monte Carlo approximation is used.
#
# For inference on sample ratios, first-order Taylor/delta covariance is
# computed from the exact unbiased moment covariance; ratio covariance itself
# is APPROXIMATE, and its estimator is not asserted unbiased.
#
# R execution status: R is unavailable in the authoring environment.
# Independent Python code was checked against exhaustive disjoint-pair
# enumeration for n=6,7,8, including marginal L2 kernels.
# Run self_test_fast_cov.R under R locally to verify this R implementation.
# For PD correction, the scaled threshold defaults to epsilon_n =
# n^(-1) * max(mean(abs(eigenvalues(n*C))), machine epsilon), which tends
# to zero provided n*C has a finite limit. Choose a larger user-supplied
# scaled threshold for more aggressive conditioning, if desired.

lcomoment_exact_cov_fast <- function(x, y, pd_floor_scaled = NULL, cov_only = FALSE) {
  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x)
  if (length(y) != n || n < 6L || any(!is.finite(x)) || any(!is.finite(y)))
    stop("Require finite equal-length vectors with n >= 6")
  if (anyDuplicated(x) || anyDuplicated(y))
    stop("Tied marginals unsupported; data must be tie-free")

  rx <- as.integer(rank(x, ties.method="first"))
  ry <- as.integer(rank(y, ties.method="first"))
  W <- cbind(1, x, y, x*x, y*y, x*y)
  P <- array(0, c(n+1L, n+1L, 6L))
  for (kk in seq_len(6L)) {
    mat <- matrix(0, n+1L, n+1L)
    mat[cbind(rx+1L, ry+1L)] <- W[, kk]
    for (ii in 2:(n+1L)) mat[ii, ] <- mat[ii, ] + mat[ii-1L, ]
    for (jj in 2:(n+1L)) mat[, jj] <- mat[, jj] + mat[, jj-1L]
    P[,,kk] <- mat
  }
  rect <- function(xlo, xhi, ylo, yhi) {
    if (xlo > xhi || ylo > yhi) return(rep(0, 6L))
    P[xhi+1L, yhi+1L, ] - P[xlo, yhi+1L, ] -
      P[xhi+1L, ylo, ] + P[xlo, ylo, ]
  }

  nm <- c("mean_X", "mean_Y", "L2_X_byY", "L2_Y_byX",
          "L3_X_byY", "L3_Y_byX", "L2_X", "L2_Y")
  orders <- c(1L, 1L, 2L, 2L, 3L, 3L, 2L, 2L)
  dimp <- length(orders)
  S0 <- numeric(dimp)
  S0[1:2] <- c(sum(x), sum(y))
  S1 <- matrix(0, n, dimp)
  S1[,1:2] <- cbind(x,y)
  T2raw <- matrix(0, dimp, dimp)
  T3raw <- matrix(0, dimp, dimp)
  ids2 <- c(3L, 4L, 7L, 8L)
  id3 <- 5:6

  for (i in 1:(n-1L)) for (j in (i+1L):n) {
    h2 <- c((x[j]-x[i])*sign(y[j]-y[i])/2,
            (y[j]-y[i])*sign(x[j]-x[i])/2,
            abs(x[j]-x[i])/2,
            abs(y[j]-y[i])/2)

    lx <- min(rx[i], rx[j]); hx <- max(rx[i], rx[j])
    ly <- min(ry[i], ry[j]); hy <- max(ry[i], ry[j])
    xlo_i <- if (rx[i] < rx[j]) i else j
    xhi_i <- if (rx[i] < rx[j]) j else i
    ylo_i <- if (ry[i] < ry[j]) i else j
    yhi_i <- if (ry[i] < ry[j]) j else i
    xb <- list(c(1L,lx-1L), c(lx+1L,hx-1L), c(hx+1L,n))
    yb <- list(c(1L,ly-1L), c(ly+1L,hy-1L), c(hy+1L,n))

    coef_a <- c(1,-2,1)
    coef_b <- c(-2*x[ylo_i]+x[yhi_i], x[ylo_i]+x[yhi_i],
                x[ylo_i]-2*x[yhi_i])
    coef_c <- c(1,-2,1)
    coef_d <- c(-2*y[xlo_i]+y[xhi_i], y[xlo_i]+y[xhi_i],
                y[xlo_i]-2*y[xhi_i])
    F <- numeric(2L)
    Q <- matrix(0, 2L, 2L)
    for (ix in 1:3) for (iy in 1:3) {
      vv <- rect(xb[[ix]][1], xb[[ix]][2],
                 yb[[iy]][1], yb[[iy]][2])
      N <- vv[1]; SX <- vv[2]; SY <- vv[3]
      SXX <- vv[4]; SYY <- vv[5]; SXY <- vv[6]
      aa <- coef_a[iy]; bb <- coef_b[iy]
      cc <- coef_c[ix]; dd <- coef_d[ix]
      F <- F + c((aa*SX+bb*N)/3, (cc*SY+dd*N)/3)
      q11 <- (aa^2*SXX + 2*aa*bb*SX + bb^2*N)/9
      q22 <- (cc^2*SYY + 2*cc*dd*SY + dd^2*N)/9
      q12 <- (aa*cc*SXY + aa*dd*SX + bb*cc*SY + bb*dd*N)/9
      Q <- Q + matrix(c(q11,q12,q12,q22),2L,2L)
    }
    S0[ids2] <- S0[ids2] + h2
    S0[id3] <- S0[id3] + F/3
    S1[i,ids2] <- S1[i,ids2] + h2
    S1[j,ids2] <- S1[j,ids2] + h2
    S1[i,id3] <- S1[i,id3] + F/2
    S1[j,id3] <- S1[j,id3] + F/2
    pair <- c(0, 0, h2[1:2], F, h2[3:4])
    T2raw <- T2raw + tcrossprod(pair)
    T3raw[id3,id3] <- T3raw[id3,id3] + Q/3
  }

  N2 <- choose(n,2L); N3 <- choose(n,3L)
  S1centered <- sweep(S1, 2L, colMeans(S1), "-")
  T1c <- crossprod(S1centered)
  S2sum <- choose(orders, 2L) * S0
  T2c <- T2raw - tcrossprod(S2sum)/N2
  T3c <- T3raw
  T3c[id3,id3] <- T3c[id3,id3] - tcrossprod(S0[id3])/N3
  num <- T1c - T2c + T3c

  C <- matrix(0, dimp, dimp)
  for (a in seq_len(dimp)) for (b in seq_len(dimp))
    C[a,b] <- num[a,b]/(choose(n,orders[a])*choose(n-orders[a],orders[b]))
  C <- (C+t(C))/2
  dimnames(C) <- list(nm, nm)
  U <- S0/choose(n,orders)
  names(U) <- nm

  # Benchmark-only option: avoid charging eigenvalue corrections and
  # nonlinear ratio inference to the covariance-computation algorithm.
  if (isTRUE(cov_only))
    return(list(n=n, moments=U, unbiased_cov=C, S1=S1))

  ev <- eigen(C, symmetric=TRUE)
  Cplus <- ev$vectors %*% (diag(pmax(ev$values,0),dimp)) %*% t(ev$vectors)
  Cplus <- (Cplus+t(Cplus))/2
  dimnames(Cplus) <- list(nm,nm)

  # Positive-definite correction is applied to S_n = n*C, as in the paper.
  # The automatic floor is sample-scale adaptive and decays as O(1/n).
  scaled_ev <- n*ev$values
  if (is.null(pd_floor_scaled)) {
    scaled_floor <- max(mean(abs(scaled_ev)), .Machine$double.eps)/n
  } else {
    scaled_floor <- as.numeric(pd_floor_scaled)
    if (length(scaled_floor) != 1L || !is.finite(scaled_floor) ||
        scaled_floor <= 0) stop("pd_floor_scaled must be a positive scalar")
  }
  # Scaling by 1/n at the end ensures the unscaled Cpd is a covariance.
  Cpd <- (ev$vectors %*%
    diag(pmax(scaled_ev, scaled_floor), dimp) %*%
    t(ev$vectors))/n
  Cpd <- (Cpd+t(Cpd))/2
  dimnames(Cpd) <- list(nm,nm)

  # Numerator is directional L-comoment; denominator is marginal L-scale.
  if (U["L2_X"] <= 0 || U["L2_Y"] <= 0)
    stop("Marginal sample L-scales must be positive")
  ind_num <- c(3L,4L,5L,6L)
  ind_den <- c(7L,8L,7L,8L)
  J <- matrix(0, 4L, dimp)
  tstat <- numeric(4L)
  for (q in 1:4) {
    aa <- ind_num[q]; bb <- ind_den[q]
    tstat[q] <- U[aa]/U[bb]
    J[q,aa] <- 1/U[bb]
    J[q,bb] <- -U[aa]/(U[bb]^2)
  }
  tnames <- c("tau2_X_byY", "tau2_Y_byX", "tau3_X_byY", "tau3_Y_byX")
  names(tstat) <- tnames
  dimnames(J) <- list(tnames, nm)
  R_unprojected <- J %*% C %*% t(J)
  R_psd <- J %*% Cplus %*% t(J)
  R_pd <- J %*% Cpd %*% t(J)
  dimnames(R_unprojected) <- list(tnames, tnames)
  dimnames(R_psd) <- list(tnames, tnames)
  dimnames(R_pd) <- list(tnames, tnames)

  list(n=n, moments=U, unbiased_cov=C, S1=S1,
       psd_corrected_cov=Cplus,
       pd_corrected_cov=Cpd,
       pd_floor_scaled=scaled_floor,
       ratio_estimates=tstat,
       ratio_cov_delta=R_unprojected,
       ratio_cov_delta_psd=R_psd,
       ratio_cov_delta_pd=R_pd,
       ratio_jacobian=J,
       minimum_eigenvalue_unbiased=min(ev$values),
       minimum_eigenvalue_pd=min(eigen(Cpd, symmetric=TRUE, only.values=TRUE)$values))
}

# Example:
# set.seed(17)
# x <- rnorm(60); y <- 0.6*x + rnorm(60)
# fit <- lcomoment_exact_cov_fast(x,y)
# fit$unbiased_cov
# fit$ratio_estimates
# fit$ratio_cov_delta_psd  # PSD first-order Taylor approximation

# fit <- lcomoment_exact_cov_fast(x, y, pd_floor_scaled = NULL)
# fit$pd_corrected_cov           # PD covariance matrix
# solve(fit$pd_corrected_cov)    # inverse exists (subject to numeric precision)
# fit$ratio_cov_delta_pd        # PD first-order ratio covariance estimate
# The scale-dependent automatic floor may be too small for near-singular data;
# set pd_floor_scaled explicitly for a stronger numerical regularization.
