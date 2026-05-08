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

  W <- fit@w
  n_genes <- nrow(W)/2
  W_pos <- W[1:n_genes,]
  W_neg <- W[(n_genes+1):(2*n_genes),]

  rownames(W_pos) <- sub("^up_", "", rownames(W_pos))
  rownames(W_neg) <- sub("^down_", "", rownames(W_neg))

  list(W = fit@w, H = fit@h, D = fit@d, W_pos = W_pos, W_neg = W_neg, fit = fit)
}

