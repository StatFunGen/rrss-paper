# Add the z-score update to a copy of the 2019 fitting loop.
make_rss2019_internal <- function(rss2019, original, output) {
  code <- paste(readLines(original), collapse = "\n")
  replace_once <- function(old, new) {
    hits <- gregexpr(old, code, fixed = TRUE)[[1]]
    stopifnot(length(hits) == 1L, hits[1] > 0)
    code <<- sub(old, new, code, fixed = TRUE)
  }
  replace_once("  attr(R, 'lambda') = lambda", "  attr(R, 'lambda') = 0")
  replace_once(
    "    s = update_each_effect_rss(R, z, s, Sigma, estimate_prior_variance,estimate_prior_method)",
    paste0("    fit_prev <- s\n",
      "    s = update_each_effect_rss(R,\n",
      "      z = drop(R %*% solve(diag(p) + R / lambda) %*% (z / lambda + colSums(s$alpha * s$mu))),\n",
      "      s, Sigma, estimate_prior_variance, estimate_prior_method)"))
  replace_once("    if((elbo[i+1]-elbo[i])<tol) {",
    "    if(norm(s$alpha - fit_prev$alpha, 'F') < tol) {")
  writeLines(code, output)
  internal <- new.env(parent = rss2019)
  sys.source(output, envir = internal)
  internal$susie_rss
}
