fit_nmf <- function(
  mat,
  k,
  n_restarts = 1,
  tol        = 1e-4,
  seed       = NULL,
  L1         = c(0., 0.)
) {
  # --- input validation ---
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)

  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_restarts, replace = FALSE)

  # --- fit ---
  fit <- RcppML::nmf(mat, k = k, seed = seeds, L1 = L1, tol = tol)

  list(W = fit@w, H = fit@h, D = fit@d)
}
