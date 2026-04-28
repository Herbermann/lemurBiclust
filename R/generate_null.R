.permute_rows_sparse <- function(mat) {
  mat <- as(mat, "RsparseMatrix")  # row-major: each row's entries are contiguous
  for (i in seq_len(nrow(mat))) {
    idx <- (mat@p[i] + 1):mat@p[i + 1]
    if (length(idx) > 1)
      mat@j[idx] <- sample(mat@j[idx])
  }
  as(mat, "CsparseMatrix")
}

.summarise_null <- function(x, probs) {
  quantile(x, probs = probs)
}


generate_null <- function(
  mat,
  fit,
  n_permutations   = 5,
  seed             = NULL,
  tol              = 1e-3,
  keep_per_factor  = FALSE,
  L1               = c(0., 0.)
) {
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)

  k       <- ncol(fit$W)
  n_genes <- nrow(mat) / 2

  probs <- c(
    seq(0,    0.90, by = 0.10),   # coarse in the bulk
    seq(0.91, 0.99, by = 0.01),   # finer in the upper tail
    seq(0.991, 0.999, by = 0.001)
  )

  set.seed(seed)

  seeds <- sample.int(.Machine$integer.max, n_permutations, replace = FALSE)

  nulls <- lapply(seq_len(n_permutations), function(i) {
    #perm_mat <- if (is(mat, "sparseMatrix")) {
    #      .permute_rows_sparse(mat)
    #    } else {
    #      t(apply(mat, 1, sample))
    #    }
    perm_mat <- t(apply(mat, 1, sample))
    #if (is(mat, "sparseMatrix"))
    #  perm_mat <- Matrix::Matrix(perm_mat, sparse = TRUE)
    
    null_fit <- RcppML::nmf(perm_mat, k = k, L1 = L1, tol = tol, seed = seeds[i])

    list(
      H     = .summarise_null(as.vector(null_fit@h), probs = probs),
      W_pos = .summarise_null(as.vector(null_fit@w[seq_len(n_genes), ]), probs = probs),
      W_neg = .summarise_null(as.vector(null_fit@w[seq(n_genes + 1, nrow(mat)), ]), probs = probs)
    )
  })

  if (keep_per_factor) {
    stop("NOT IMPLEMENTED")

  } else {
    list(
      probs = probs,
      H     = rowMeans(do.call(cbind, lapply(nulls, `[[`, "H"))),
      W_pos = rowMeans(do.call(cbind, lapply(nulls, `[[`, "W_pos"))),
      W_neg = rowMeans(do.call(cbind, lapply(nulls, `[[`, "W_neg")))
    )
  }
}