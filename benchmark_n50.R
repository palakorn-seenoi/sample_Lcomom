# ================================================================
# L-comoment exact unbiased covariance: direct vs. fast algorithm
# Sample size n=50; Direct algorithm timed exactly once.
#
# Generated for the paper:
# Exact Finite-Sample Covariances of Second- and Third-Order L-Comoments
#
# HOW TO RUN from R / RStudio, after extracting both .R files:
#   source("D:/justwork/Lcomoment_Benchmark_n50/benchmark_n50.R", chdir=TRUE)
# (Adjust the path to the actual unzipped folder.)
# No external R packages are required. Outputs are saved in the
# directory containing the script if source(..., chdir=TRUE).
# ================================================================

N <- 50L
SEED <- 20261017L
RHO <- 0.60
FAST_REPS <- 3L
DIRECT_REPS <- 1L
FAST_BATCH <- 3L
DIRECT_BATCH <- 1L
RUN_SELF_TEST <- TRUE

if (!file.exists("lcomoment_fast_exact_cov.R")) {
  stop("Cannot find lcomoment_fast_exact_cov.R in the working directory. ",
       "Unzip both R files to one folder, then call source(..., chdir=TRUE).")
}
source("lcomoment_fast_exact_cov.R")
fast_cov_only <- function(x, y)
  lcomoment_exact_cov_fast(x, y, cov_only=TRUE)

# ---------------------------------------------------------------
# Literal disjoint-subset reference using all 1-, 2-, and 3-subsets.
# Internally vectorized over the second subset to limit R loop
# overhead. Still enumerates the disjoint combinations (O(n^6)
# work for third-order/third-order covariance).
# Output is the FULL 8x8 covariance, matching the fast function.
# ---------------------------------------------------------------
lcomoment_exact_cov_direct <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  n <- length(x)
  if (n != length(y) || n < 6L ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Need finite x,y, equal length, n >= 6.")
  }
  if (anyDuplicated(x) || anyDuplicated(y))
    stop("No ties allowed in x or y.")

  # Columns are index subsets; each column has m indices.
  idx <- list(matrix(seq_len(n), nrow=1L),
              combn(n, 2L), combn(n, 3L))

  i2 <- idx[[2L]]
  xa <- x[i2[1L, ]]; xb <- x[i2[2L, ]]
  ya <- y[i2[1L, ]]; yb <- y[i2[2L, ]]
  H2 <- cbind(
    (xb-xa)*sign(yb-ya)/2,
    (yb-ya)*sign(xb-xa)/2,
    abs(xb-xa)/2,
    abs(yb-ya)/2
  )

  i3 <- idx[[3L]]
  third_xy <- vapply(seq_len(ncol(i3)), function(k) {
    ids <- i3[,k]
    z <- x[ids][order(y[ids])]
    (z[1L]-2*z[2L]+z[3L])/3
  }, numeric(1))
  third_yx <- vapply(seq_len(ncol(i3)), function(k) {
    ids <- i3[,k]
    z <- y[ids][order(x[ids])]
    (z[1L]-2*z[2L]+z[3L])/3
  }, numeric(1))

  H <- list(cbind(x, y), H2, cbind(third_xy, third_yx))
  # Moment order in fast routine:
  # mean_X,mean_Y,L2_X_byY,L2_Y_byX,L3_X_byY,L3_Y_byX,L2_X,L2_Y
  slots <- list(c(1L,2L), c(3L,4L,7L,8L), c(5L,6L))
  nm <- c("mean_X", "mean_Y", "L2_X_byY", "L2_Y_byX",
          "L3_X_byY", "L3_Y_byX", "L2_X", "L2_Y")

  # contain[[s]][j,i] is TRUE if subset j of order s contains i.
  contain <- lapply(idx, function(I) {
    ns <- ncol(I)
    A <- matrix(FALSE, nrow=ns, ncol=n)
    for (k in seq_len(nrow(I))) {
      A[cbind(seq_len(ns), I[k,])] <- TRUE
    }
    A
  })

  C <- matrix(0, 8L, 8L, dimnames=list(nm,nm))
  moments <- numeric(8L)
  names(moments) <- nm
  for (m in 1:3) moments[slots[[m]]] <- colMeans(H[[m]])

  for (m in 1:3) for (s in m:3) {
    Hm <- H[[m]]
    Hs <- H[[s]]
    Im <- idx[[m]]
    ns_m <- ncol(Im)
    ns_s <- ncol(idx[[s]])
    weighted_sum <- matrix(0, nrow=ncol(Hm), ncol=ncol(Hs))

    # Explicitly enumerate the first index subset I. For each I,
    # scan all J of order s and retain only J disjoint from I.
    for (j in seq_len(ns_m)) {
      shared <- rowSums(contain[[s]][, Im[,j], drop=FALSE])
      admissible <- shared == 0L
      sum_J <- colSums(Hs[admissible, , drop=FALSE])
      weighted_sum <- weighted_sum + outer(Hm[j,], sum_J)
    }

    n_disjoint_pairs <- ns_m * choose(n-m, s)
    D <- weighted_sum / n_disjoint_pairs
    # Exact unbiased covariance, not asymptotic approximation.
    block <- outer(colMeans(Hm), colMeans(Hs)) - D
    C[slots[[m]], slots[[s]]] <- block
    if (s != m) C[slots[[s]], slots[[m]]] <- t(block)
  }
  C <- (C + t(C))/2
  list(n=n, moments=moments, unbiased_cov=C)
}

compare_covariances <- function(fast, direct, n, tol=1e-9) {
  A <- fast$unbiased_cov
  B <- direct$unbiased_cov
  dif <- abs(A-B)
  max_abs_diff <- max(dif)
  scale <- max(1, max(abs(A)), max(abs(B)))
  max_moment_diff <- max(abs(fast$moments-direct$moments))
  cat(sprintf("n=%d | max matrix error=%.5g | max moment error=%.5g\n",
              n, max_abs_diff, max_moment_diff))
  if (!is.finite(max_abs_diff) ||
      max_abs_diff > tol*scale || max_moment_diff > tol*scale) {
    stop("Direct and fast algorithms disagree. Please send output for debugging.")
  }
  invisible(list(max_abs_diff=max_abs_diff,
                 max_moment_diff=max_moment_diff))
}

# Reproducibility check at a very small n before benchmarking.
if (RUN_SELF_TEST) {
  cat("=== Preliminary correctness test (n=7, n=8) ===\n")
  for (ncheck in c(7L,8L)) {
    set.seed(SEED+ncheck)
    x0 <- rnorm(ncheck)
    y0 <- RHO*x0 + sqrt(1-RHO^2)*rnorm(ncheck)
    fit0 <- fast_cov_only(x0, y0)
    ref0 <- lcomoment_exact_cov_direct(x0, y0)
    compare_covariances(fit0, ref0, ncheck)
  }
}

# Fixed simulated dataset: SAME observations for both methods.
set.seed(SEED)
x <- rnorm(N)
y <- RHO*x + sqrt(1-RHO^2)*rnorm(N)
cat("\n=== Benchmark setup ===\n")
cat(sprintf("N=%d, seed=%d, underlying Gaussian rho=%.2f\n", N, SEED, RHO))
cat(sprintf("Number of disjoint 3+3 subset pairs: %s\n",
            format(choose(N,3)*choose(N-3,3), big.mark=",")))
cat("Runtime measures FULL 8x8 unbiased covariance matrices.\n")

# Fast warm-up is excluded from measured times.
cat("\n=== Untimed Fast warm-up ===\n")
fast_fit <- fast_cov_only(x, y)

benchmark <- function(label, fun, x, y, reps, batch) {
  measurements <- numeric(reps)
  last_answer <- NULL
  for (r in seq_len(reps)) {
    t <- system.time({
      for (b in seq_len(batch)) last_answer <- fun(x, y)
    })[["elapsed"]]
    measurements[r] <- t / batch
    cat(sprintf("%-11s repetition %d/%d: %.6f sec per full 8x8 matrix\n",
                label, r, reps, measurements[r]))
  }
  list(times=measurements, answer=last_answer)
}

cat("\n=== Timings: Fast 3 repetitions, Direct 1 repetition ===\n")
fast_result <- benchmark("Fast", fast_cov_only,
                         x, y, FAST_REPS, FAST_BATCH)

# IMPORTANT: There is no untimed Direct computation. Reuse the
# result of this sole Direct call for the accuracy comparison.
cat("\nStarting the one Direct computation (n=50).\n")
direct_result <- benchmark("Direct", lcomoment_exact_cov_direct,
                           x, y, DIRECT_REPS, DIRECT_BATCH)
fast_t <- fast_result$times
direct_t <- direct_result$times
fast_fit <- fast_result$answer
direct_fit <- direct_result$answer

cat("\n=== Full-matrix numerical comparison ===\n")
check <- compare_covariances(fast_fit, direct_fit, N)
fast_median <- median(fast_t)
direct_median <- median(direct_t)
speedup <- if (fast_median > 0) direct_median / fast_median else NA_real_

entries <- data.frame(
  quantity=c("Var(L3_X_byY)", "Cov(L3_X_byY,L3_Y_byX)", "Var(L3_Y_byX)"),
  direct=c(direct_fit$unbiased_cov["L3_X_byY","L3_X_byY"],
           direct_fit$unbiased_cov["L3_X_byY","L3_Y_byX"],
           direct_fit$unbiased_cov["L3_Y_byX","L3_Y_byX"]),
  fast=c(fast_fit$unbiased_cov["L3_X_byY","L3_X_byY"],
         fast_fit$unbiased_cov["L3_X_byY","L3_Y_byX"],
         fast_fit$unbiased_cov["L3_Y_byX","L3_Y_byX"]),
  stringsAsFactors=FALSE
)
entries$abs_difference <- abs(entries$direct-entries$fast)

summary_row <- data.frame(
  n=N, seed=SEED, rho=RHO,
  number_disjoint_triple_pairs=choose(N,3)*choose(N-3,3),
  direct_median_sec=direct_median,
  fast_median_sec=fast_median,
  speedup_direct_over_fast=speedup,
  max_absolute_matrix_error=check$max_abs_diff,
  max_absolute_moment_error=check$max_moment_diff,
  repetitions=DIRECT_REPS,
  fast_repetitions=FAST_REPS,
  direct_repetitions=DIRECT_REPS,
  stringsAsFactors=FALSE
)
raw_timing <- rbind(
  data.frame(method="Direct", repetition=seq_len(DIRECT_REPS),
             seconds_per_full_matrix=direct_t),
  data.frame(method="Fast", repetition=seq_len(FAST_REPS),
             seconds_per_full_matrix=fast_t)
)

prefix <- paste0("lcomoment_n",N,"_")
write.csv(summary_row, paste0(prefix,"timing_summary.csv"), row.names=FALSE)
write.csv(raw_timing, paste0(prefix,"timing_repetitions.csv"), row.names=FALSE)
write.csv(entries, paste0(prefix,"covariance_examples.csv"), row.names=FALSE)
write.csv(direct_fit$unbiased_cov,
          paste0(prefix,"direct_full_cov.csv"), row.names=TRUE)
write.csv(fast_fit$unbiased_cov,
          paste0(prefix,"fast_full_cov.csv"), row.names=TRUE)

cat("\n=== Summary (actual measured values) ===\n")
print(summary_row, row.names=FALSE, digits=6)
cat("\nThird-order covariance entries:\n")
print(entries, row.names=FALSE, digits=12)
cat("\nFiles saved in: ",normalizePath(getwd()),"\n",sep="")
cat("Files: timing_summary, timing_repetitions, covariance_examples, ",
    "direct_full_cov, fast_full_cov (CSV).\n",sep="")
cat("NOTE: n=50 Direct time is ONE measurement, not a 3-run median.\n")
