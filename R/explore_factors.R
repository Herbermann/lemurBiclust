genespace_stability <- function(fit1, fit2) {
  # weight W columns by their factor importance
  W1 <- fit1@w %*% diag(fit1@d)
  W2 <- fit2@w %*% diag(fit2@d)
  
  # remove dead factors
  W1 <- W1[, colSums(W1) > 0, drop = FALSE]
  W2 <- W2[, colSums(W2) > 0, drop = FALSE]
  
  # now QR on the weighted W
  Q1 <- qr.Q(qr(W1))
  Q2 <- qr.Q(qr(W2))
  
  decomp <- RcppML::svd(t(Q1) %*% Q2)
  theta2 <- acos(decomp@d)^2
  sqrt(theta2)/(length(theta2)* 0.5 * pi)
}



explore_factors <- function(
    mat,
    ks,
    n_restarts      = 2,
    test_fraction   = 0.,
    L1              = c(0., 0.),
    tol             = 1e-4,
    show_plot       = FALSE
) {

  # --- validate input ---
  if (!is.numeric(ks) || any(ks < 1) || any(ks != floor(ks)))
    stop("'ks' must be a vector of positive integers")
  if (any(ks >= min(dim(mat))))
    stop("all 'ks' must be smaller than both dimensions of 'mat'")
  if (!is.numeric(n_restarts) || length(n_restarts) != 1 || n_restarts < 1)
    stop("'n_restarts' must be a single positive integer")
  if (!is.numeric(test_fraction) || test_fraction < 0 || test_fraction >= 1)
    stop("'test_fraction' must be a number in [0, 1)")

  # -- grab a progressbar ---
  p <- progressr::progressor(steps = length(ks) * n_restarts)


  # --- Preps for validation set --
  use_cv <- test_fraction > 0

  if (use_cv) {
    # eligible pool computed once; mask resampled per restart
    nz_idx <- which(mat != 0, arr.ind = TRUE)
    n_test <- round(nrow(nz_idx) * test_fraction)
  }

  # --- main logic ---
  results <- lapply(ks, function(k) {

    fits_and_errors <- lapply(seq_len(n_restarts), function(i) {

      if (use_cv) {
        # fresh mask per restart
        test_coords  <- nz_idx[sample.int(nrow(nz_idx), n_test), , drop = FALSE]
        test_mask    <- Matrix::sparseMatrix(
          i    = test_coords[, 1],
          j    = test_coords[, 2],
          x    = 1,
          dims = dim(mat)
        )
        train_mat <- mat
        train_mat[test_coords] <- 0
        if (inherits(mat, "sparseMatrix"))
          train_mat <- methods::as(train_mat, "dgCMatrix")
      } else {
        train_mat <- mat
      }

      fit <- RcppML::nmf(train_mat, k = k, seed = i, L1 = L1, tol = tol)

      err <- if (use_cv) {
        RcppML::evaluate(fit, data = mat, mask = test_mask, missing_only = TRUE)
      } else {
        fit@misc[["loss"]]
      }

      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))

      list(fit = fit, error = err)
    })

    fits   <- lapply(fits_and_errors, `[[`, "fit")
    errors <- sapply(fits_and_errors, `[[`, "error")
    best_idx <- which.min(errors)

    # --- stability -----------------------------------------------------------
    if (length(fits) < 2) {
      sims <- NA_real_
    } else {
      sims <- combn(length(fits), 2, function(idx) {
        tryCatch({
          if (ncol(fits[[idx[1]]]@w) == 0 || ncol(fits[[idx[2]]]@w) == 0) return(NA_real_)
          mean(genespace_stability(fits[[idx[1]]], fits[[idx[2]]]))
        }, error = function(e) NA_real_)
      }, simplify = TRUE)
    }

    list(
      k         = k,
      error     = mean(errors),
      stability = max(sims, na.rm = TRUE),
      best_fit  = fits[[best_idx]]
    )
  })

  error_df <- data.frame(k = ks, error     = sapply(results, `[[`, "error"))
  stab_df  <- data.frame(k = ks, stability = sapply(results, `[[`, "stability"))

  p_error <- ggplot(error_df, aes(x = k, y = error)) +
    geom_line() + geom_point() +
    labs(x = "k", y = if (use_cv) "test MSE" else "reconstruction error") +
    theme_classic()

  p_stab <- ggplot(stab_df, aes(x = k, y = stability)) +
    geom_line() + geom_point() +
    labs(x = "k", y = "mean factor stability") +
    theme_classic()

  structure(
    list(
      results  = results,
      best_fits = lapply(results, function(r) r$best_fit),  # see below
      summary  = data.frame(
        k         = ks,
        error     = sapply(results, `[[`, "error"),
        stability = sapply(results, `[[`, "stability")
      ),
      plot     = patchwork::wrap_plots(p_error, p_stab),
      params   = list(L1 = L1, tol = tol, n_restarts = n_restarts,
                      test_fraction = test_fraction)
    ),
    class = "factor_exploration"
  )
}

explore_factors_progress <- function(..., handler = "cli") {
  progressr::handlers(handler)
  progressr::with_progress(explore_factors(...))
}

print.factor_exploration <- function(x, ...) {
  cat("Factor exploration over k =", paste(x$summary$k, collapse = ", "), "\n")
  print(x$summary)
  invisible(x)
}

plot.factor_exploration <- function(x, ...) {
  print(x$plot)
  invisible(x)
}