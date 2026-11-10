args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1])
out <- normalizePath(args[2])
commit <- args[3]
suppressPackageStartupMessages({library(Matrix); library(matrixStats)})
original <- file.path(out, "historical_source")
dir.create(original, showWarnings = FALSE)
utils::untar(file.path(out, "susieR.tar.gz"), exdir = original)
source_dir <- file.path(original, paste0("susieR-", commit))
stopifnot(read.dcf(file.path(source_dir, "DESCRIPTION"))[1, "Version"] == "0.7.1.0487")
rss2019 <- new.env(parent = globalenv())
for (f in list.files(file.path(source_dir, "R"), pattern = "[.]R$", full.names = TRUE))
  sys.source(f, envir = rss2019)

rss2019_in_new_bottle <- function(z, R, lambda, L = 3,
                           max_iter = 100, tol = 1e-5) {
  p <- length(z)
  fit <- list(alpha = matrix(0, L, p), mu = matrix(0, L, p))
  for (i in seq_len(max_iter)) {
    fit_prev <- fit
    fit <- rss2019$susie_rss(
      z = drop(R %*% solve(diag(p) + R / lambda) %*% (z / lambda + colSums(fit$alpha * fit$mu))),
      R = R, L = L, s_init = fit,
      prior_variance = 1e20,
      estimate_residual_variance = FALSE, estimate_prior_variance = FALSE,
      max_iter = 1, coverage = NULL)
    if (norm(fit$alpha - fit_prev$alpha, "F") < tol) break
  }
  fit
}

source(file.path(root, "rss2019_internal.R"))
rss_internal <- make_rss2019_internal(rss2019,
  file.path(source_dir, "R", "susie_rss.R"),
  file.path(out, "susie_rss_internal.R"))
cases <- jsonlite::read_json(file.path(out, "inputs.json"), simplifyVector = FALSE)
fits <- lapply(cases, function(case) {
  z <- unlist(case$z)
  R <- do.call(rbind, lapply(case$R, unlist))
  lambda <- case$a
  fit <- rss2019_in_new_bottle(z, R, lambda)
  internal <- rss_internal(z, R, L = 3, lambda = lambda,
    prior_variance = 1e20,
    estimate_prior_variance = FALSE, estimate_residual_variance = FALSE,
    tol = 1e-5, coverage = NULL,
    s_init = list(alpha = matrix(0, 3, ncol(R)), mu = matrix(0, 3, ncol(R))))
  stopifnot(max(abs(fit$alpha - internal$alpha)) < 1e-10)
  list(case = case$case, alpha = fit$alpha, internal_alpha = internal$alpha)
})
jsonlite::write_json(fits, file.path(out, "rss_fits.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
