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
      k           = k,
      errors      = errors,
      stabilities = sims,
      error       = mean(errors),
      stability   = max(sims, na.rm = TRUE),
      best_fit    = fits[[best_idx]]
    )
  })

  # --- prepare data frames for plotting.
  error_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      k = x$k,
      error = x$errors
    )
  }))

  error_summary <- error_df |>
    dplyr::summarise(
      mean = mean(error),
      sd = sd(error),
      .by = k
    )

  stab_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      k = rep(x$k, length(x$stabilities)),
      stability = x$stabilities
    )
  }))

  stab_summary <- stab_df |>
    dplyr::summarise(
      mean = mean(stability, na.rm = TRUE),
      sd   = sd(stability, na.rm = TRUE),
      .by = k
  )

  # --- error plot
  p_error <-
    ggplot(error_df, aes(k, error)) +
    geom_errorbar(
      data = error_summary,
      aes(
        x = k,
        ymin = mean - sd,
        ymax = mean + sd
      ),
      inherit.aes = FALSE,
      width = 0.15,
      linewidth = 0.5
    ) + 
    geom_ribbon(
      data = error_summary,
      aes(
        x = k,
        ymin = mean - sd,
        ymax = mean + sd
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    )+
    geom_line(data = error_summary, aes(y = mean)) +
    geom_point(data = error_summary, aes(y = mean)) +
    scale_x_continuous(breaks = ks) + 
    labs(
      x = "Factor number",
      y = if (use_cv) "Reconstruction error" else "Reconstruction error"
    ) +
    theme_classic(base_size = 14)


  # --- stability plot
  p_stab <-
  ggplot(stab_df, aes(k, stability)) +
  geom_errorbar(
    data = stab_summary,
    aes(
      x = k,
      ymin = mean - sd,
      ymax = mean + sd
    ),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 0.5
  ) + 
  geom_ribbon(
    data = stab_summary,
    aes(
      x = k,
      ymin = mean - sd,
      ymax = mean + sd
    ),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = stab_summary,
    aes(y = mean),
    linewidth = 1
  ) +
  geom_point(
    data = stab_summary,
    aes(y = mean),
    size = 2
  ) +
  scale_x_continuous(breaks = ks) +
  labs(
    x = "Factor number",
    y = "Factor stability"
  ) +
  theme_classic(base_size = 14)
  
  # --- output structure.
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



############

genespace_stability <- function(fit1, fit2) {
  # weight W columns by factor importance
  W1 <- fit1@w %*% diag(fit1@d)
  W2 <- fit2@w %*% diag(fit2@d)

  # remove dead factors
  W1 <- W1[, colSums(W1) > 0, drop = FALSE]
  W2 <- W2[, colSums(W2) > 0, drop = FALSE]

  # principal angles between factor subspaces
  Q1 <- qr.Q(qr(W1))
  Q2 <- qr.Q(qr(W2))

  decomp <- RcppML::svd(t(Q1) %*% Q2)
  theta2 <- acos(decomp@d)^2

  sqrt(theta2) / (length(theta2) * 0.5 * pi)
}

normalize_l1 <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 1) x <- c(x, x)
  if (length(x) != 2) stop("Each L1 setting must have length 1 or 2.")
  x
}

format_l1 <- function(x) {
  sprintf("c(%g, %g)", x[1], x[2])
}
# helper: sparsity of a fitted factorization
factor_sparsity_summary <- function(fit, zero_tol = 1e-8) {
  W <- fit@w
  H <- fit@h

  data.frame(
    sparsity_w = mean(abs(W) <= zero_tol, na.rm = TRUE),
    sparsity_h = mean(abs(H) <= zero_tol, na.rm = TRUE)
  )
}


pareto_front <- function(df, x, y) {
  # x minimized, y maximized

  ord <- order(df[[x]], -df[[y]])

  best_y <- -Inf
  keep <- logical(nrow(df))

  for (i in ord) {
    yi <- df[[y]][i]

    if (yi > best_y) {
      keep[i] <- TRUE
      best_y <- yi
    }
  }

  keep
}


explore_sparsity <- function(
    mat,
    k,
    l1_grid,
    n_restarts   = 2,
    test_fraction = 0,
    tol           = 1e-4,
    show_plot     = FALSE
) {

  # --- validate input ---
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
    stop("'k' must be a single positive integer")
  }

  if (k >= min(dim(mat))) {
    stop("'k' must be smaller than both dimensions of 'mat'")
  }

  if (!is.numeric(n_restarts) || length(n_restarts) != 1 || n_restarts < 1) {
    stop("'n_restarts' must be a single positive integer")
  }

  if (!is.numeric(test_fraction) || test_fraction < 0 || test_fraction >= 1) {
    stop("'test_fraction' must be a number in [0, 1)")
  }

  if (is.atomic(l1_grid) && !is.list(l1_grid)) {
    l1_grid <- as.list(l1_grid)
  }

  if (!is.list(l1_grid) || length(l1_grid) == 0) {
    stop("'l1_grid' must be a non-empty list or numeric vector")
  }

  l1_grid <- lapply(l1_grid, normalize_l1)
  l1_labels <- vapply(l1_grid, format_l1, character(1))

  use_cv <- test_fraction > 0

  # --- CV masks ---
  cv_masks <- NULL

  if (use_cv) {

    nz_idx <- which(mat != 0, arr.ind = TRUE)

    if (nrow(nz_idx) == 0) {
      stop("'mat' has no non-zero entries for CV masking")
    }

    n_test <- max(1, round(nrow(nz_idx) * test_fraction))

    cv_masks <- lapply(seq_len(n_restarts), function(i) {

      test_coords <- nz_idx[
        sample.int(nrow(nz_idx), n_test),
        ,
        drop = FALSE
      ]

      test_mask <- Matrix::sparseMatrix(
        i = test_coords[, 1],
        j = test_coords[, 2],
        x = 1,
        dims = dim(mat)
      )

      list(coords = test_coords, mask = test_mask)
    })
  }

  p <- progressr::progressor(
    steps = length(l1_grid) * n_restarts
  )

  # --- run fits ---
  results <- lapply(seq_along(l1_grid), function(li) {

    L1 <- l1_grid[[li]]

    fits_and_errors <- lapply(seq_len(n_restarts), function(i) {

      if (use_cv) {

        test_coords <- cv_masks[[i]]$coords
        test_mask <- cv_masks[[i]]$mask

        train_mat <- mat
        train_mat[test_coords] <- 0

        if (inherits(train_mat, "sparseMatrix")) {
          train_mat <- methods::as(train_mat, "dgCMatrix")
        }

      } else {

        train_mat <- mat
        test_mask <- NULL
      }

      fit <- RcppML::nmf(
        train_mat,
        k = k,
        seed = i,
        L1 = L1,
        tol = tol
      )

      err <- if (use_cv) {

        RcppML::evaluate(
          fit,
          data = mat,
          mask = test_mask,
          missing_only = TRUE
        )

      } else {

        fit@misc[["loss"]]
      }

      p(
        message = sprintf(
          "L1=%s, restart %d/%d",
          l1_labels[li],
          i,
          n_restarts
        )
      )

      list(
        fit = fit,
        error = err
      )
    })

    fits <- lapply(fits_and_errors, `[[`, "fit")
    errors <- sapply(fits_and_errors, `[[`, "error")

    best_idx <- which.min(errors)

    sps <- do.call(
      rbind,
      lapply(fits, factor_sparsity_summary)
    )

    mean_sparsity_w <- mean(sps$sparsity_w, na.rm = TRUE)
    mean_sparsity_h <- mean(sps$sparsity_h, na.rm = TRUE)

    list(
      L1 = L1,
      label = l1_labels[li],
      error = mean(errors, na.rm = TRUE),
      sparsity_w = mean_sparsity_w,
      sparsity_h = mean_sparsity_h,
      best_fit = fits[[best_idx]],
      fits = fits,
      errors = errors
    )
  })

  # --- summary ---
  summary_df <- data.frame(
    L1_label = vapply(results, `[[`, character(1), "label"),
    L1_w = vapply(results, function(x) x$L1[1], numeric(1)),
    L1_h = vapply(results, function(x) x$L1[2], numeric(1)),
    error = vapply(results, `[[`, numeric(1), "error"),
    sparsity_w = vapply(results, `[[`, numeric(1), "sparsity_w"),
    sparsity_h = vapply(results, `[[`, numeric(1), "sparsity_h"),
    stringsAsFactors = FALSE
  )

  # --- Pareto fronts ---
  summary_df$pareto_w <- pareto_front(
    summary_df,
    x = "error",
    y = "sparsity_w"
  )

  summary_df$pareto_h <- pareto_front(
    summary_df,
    x = "error",
    y = "sparsity_h"
  )

  # --- plots ---
  p_error <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = L1_label,
      y = error,
      group = 1
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "L1 setting",
      y = if (use_cv) {
        "held-out loss"
      } else {
        "reconstruction loss"
      }
    ) +
    ggplot2::theme_classic()

  p_sw <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = L1_label,
      y = sparsity_w,
      group = 1
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "L1 setting",
      y = "mean sparsity of W"
    ) +
    ggplot2::theme_classic()

  p_sh <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = L1_label,
      y = sparsity_h,
      group = 1
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "L1 setting",
      y = "mean sparsity of H"
    ) +
    ggplot2::theme_classic()

  # --- Pareto plots ---
  p_pareto_w <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = error,
      y = sparsity_w,
      label = L1_label
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_text(
      nudge_y = 0.01,
      size = 3
    ) +
    ggplot2::geom_point(
      data = summary_df[summary_df$pareto_w, ],
      size = 3
    ) +
    ggplot2::geom_path(
      data = summary_df[
        order(summary_df$error) &
          summary_df$pareto_w,
      ]
    ) +
    ggplot2::labs(
      x = "reconstruction error",
      y = "sparsity of W",
      title = "Pareto frontier: W"
    ) +
    ggplot2::theme_classic()

  p_pareto_h <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = error,
      y = sparsity_h,
      label = L1_label
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_text(
      nudge_y = 0.01,
      size = 3
    ) +
    ggplot2::geom_point(
      data = summary_df[summary_df$pareto_h, ],
      size = 3
    ) +
    ggplot2::geom_path(
      data = summary_df[
        order(summary_df$error) &
          summary_df$pareto_h,
      ]
    ) +
    ggplot2::labs(
      x = "reconstruction error",
      y = "sparsity of H",
      title = "Pareto frontier: H"
    ) +
    ggplot2::theme_classic()

  out <- list(
    results = results,
    best_fits = lapply(results, `[[`, "best_fit"),
    summary = summary_df,
    plot = patchwork::wrap_plots(
      p_error,
      p_sw,
      p_sh,
      p_pareto_w,
      p_pareto_h
    ),
    params = list(
      k = k,
      l1_grid = l1_grid,
      n_restarts = n_restarts,
      test_fraction = test_fraction,
      tol = tol
    )
  )

  class(out) <- "sparsity_exploration"

  out
}

explore_sparsity_progress <- function(..., handler = "cli") {
  progressr::handlers(handler)
  progressr::with_progress(explore_sparsity(...))
}

print.sparsity_exploration <- function(x, ...) {
  cat("Sparsity exploration at fixed k =", x$params$k, "\n")
  print(x$summary)
  invisible(x)
}

plot.sparsity_exploration <- function(x, ...) {
  print(x$plot)
  invisible(x)
}