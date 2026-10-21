# ================================================================
# Monte Carlo validation of exact unbiased L-comoment covariance
# 8-dimensional vector; no external packages required.
# Run: source("D:/justwork/Lcomoment_MC_validation/monte_carlo_lcomoment_cov.R", chdir=TRUE)
# ================================================================
N_VALUES <- c(20L, 30L, 50L)
B <- 1000L                  # Pilot: 1000; optionally increase to 3000 later
SEED <- 20261018L
RHO <- 0.60
PRINT_EVERY <- 100L
OUTPUT_DIR <- getwd()

source("lcomoment_fast_exact_cov.R")
stopifnot(B >= 30L, all(N_VALUES >= 6L), abs(RHO) < 1)

nm <- c("mean_X", "mean_Y", "L2_X_byY", "L2_Y_byX",
        "L3_X_byY", "L3_Y_byX", "L2_X", "L2_Y")
p <- length(nm)

# Per-cell Monte Carlo uncertainty for difference mean(C_hat) - cov(U).
# The per-rep pseudo difference uses empirical centering and produces
# a first-order valid Monte Carlo SE (no assertion of independence
# between C_hat and outer product U*U).
run_n <- function(n, B) {
  cat(sprintf("\n=== Monte Carlo n=%d, B=%d, rho=%.2f ===\n", n, B, RHO))
  U <- matrix(NA_real_, nrow=B, ncol=p, dimnames=list(NULL,nm))
  C_arr <- array(NA_real_, dim=c(p,p,B), dimnames=list(nm,nm,NULL))
  min_ev <- numeric(B)
  t0 <- proc.time()[["elapsed"]]
  for (b in seq_len(B)) {
    x <- rnorm(n)
    y <- RHO*x + sqrt(1-RHO^2)*rnorm(n)
    fit <- lcomoment_exact_cov_fast(x, y, cov_only=TRUE)
    U[b, ] <- fit$moments[nm]
    C_arr[,,b] <- fit$unbiased_cov[nm,nm]
    min_ev[b] <- min(eigen(fit$unbiased_cov, symmetric=TRUE,
                           only.values=TRUE)$values)
    if (b %% PRINT_EVERY == 0L || b == B) {
      cat(sprintf(" n=%d: %d/%d complete (%.1f seconds)\n",
                  n, b, B, proc.time()[["elapsed"]]-t0))
      flush.console()
    }
  }
  mean_C <- apply(C_arr, c(1L,2L), mean)
  reference_C <- stats::cov(U)  # empirical covariance of the sample L-comoment vector
  delta <- mean_C-reference_C
  centered_U <- sweep(U, 2L, colMeans(U), "-")
  # Approximate influence values: C_b - (U_b - E U)(U_b - E U)^T
  # Using sample mean in place of E U is first-order equivalent.
  infl <- matrix(NA_real_, B, p*p)
  for (b in seq_len(B))
    infl[b,] <- as.vector(C_arr[,,b] - tcrossprod(centered_U[b,]))
  se_delta <- matrix(apply(infl, 2L, sd)/sqrt(B), p, p,
                     dimnames=list(nm,nm))
  z <- delta/se_delta
  z[se_delta == 0] <- NA_real_
  # Diagonal sample SD of C_hat across simulations; estimation variability
  sd_C <- matrix(apply(matrix(aperm(C_arr,c(3,1,2)), nrow=B), 2, sd),p,p,
                 dimnames=list(nm,nm))
  ix <- which(upper.tri(mean_C, diag=TRUE), arr.ind=TRUE)
  detail <- data.frame(n=n, row_name=nm[ix[,1]], col_name=nm[ix[,2]],
                       mean_estimated_cov=mean_C[ix],
                       mc_reference_cov=reference_C[ix],
                       difference=delta[ix], mc_se_difference=se_delta[ix],
                       z_difference=z[ix], sd_estimated_cov=sd_C[ix],
                       stringsAsFactors=FALSE)
  neg <- min_ev < -1e-10
  summary <- data.frame(n=n, B=B, rho=RHO, seed=SEED,
                        max_abs_difference=max(abs(delta)),
                        max_abs_z=max(abs(z),na.rm=TRUE),
                        rms_difference=sqrt(mean(delta^2)),
                        relative_frobenius_difference=sqrt(sum(delta^2))/sqrt(sum(reference_C^2)),
                        non_psd_count=sum(neg), non_psd_rate=mean(neg),
                        min_eigenvalue_median=median(min_ev),
                        min_eigenvalue_min=min(min_ev),
                        elapsed_sec=proc.time()[["elapsed"]]-t0)
  # Save a concise covariance matrix comparison by n.
  write.csv(mean_C, file.path(OUTPUT_DIR,sprintf("n%d_mean_unbiased_cov.csv",n)))
  write.csv(reference_C,file.path(OUTPUT_DIR,sprintf("n%d_mc_reference_cov.csv",n)))
  write.csv(detail,file.path(OUTPUT_DIR,sprintf("n%d_entrywise_comparison.csv",n)),row.names=FALSE)
  write.csv(data.frame(min_eigenvalue=min_ev, non_psd=neg),
            file.path(OUTPUT_DIR,sprintf("n%d_eigenvalues.csv",n)),row.names=FALSE)
  write.csv(summary,file.path(OUTPUT_DIR,sprintf("n%d_summary.csv",n)),row.names=FALSE)
  print(summary, row.names=FALSE, digits=5)
  list(summary=summary,detail=detail)
}

set.seed(SEED)
all_results <- lapply(N_VALUES, run_n, B=B)
summary_all <- do.call(rbind,lapply(all_results,`[[`,"summary"))
detail_all <- do.call(rbind,lapply(all_results,`[[`,"detail"))
write.csv(summary_all,file.path(OUTPUT_DIR,"mc_summary_all.csv"),row.names=FALSE)
write.csv(detail_all,file.path(OUTPUT_DIR,"mc_entrywise_all.csv"),row.names=FALSE)
cat("\n=== FINAL SUMMARY ===\n")
print(summary_all,row.names=FALSE,digits=6)
cat("\nSaved results to:",normalizePath(OUTPUT_DIR),"\n")
cat("Please upload mc_summary_all.csv and mc_entrywise_all.csv first.\n")
