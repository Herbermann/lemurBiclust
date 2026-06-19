library(RcppML)
library(clue)
library(withr)

`%||%` <- function(x, y) if (!is.null(x)) x else y

.cosine_sim <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  nx <- sqrt(sum(x^2))
  ny <- sqrt(sum(y^2))

  if (nx == 0 || ny == 0) return(0)

  sum(x * y) / (nx * ny)
}

.get_wh <- function(fit) {
  w <- fit$w %||% fit$W
  h <- fit$h %||% fit$H

  if (is.null(w) || is.null(h)) {
    stop("Could not find factor matrices in the fit object. Expected $w/$h or $W/$H.")
  }

  list(
    w = as.matrix(w),
    h = as.matrix(h)
  )
}

.compose_nmf_input <- function(X, ...) {
  if (inherits(X, "Matrix") || is.matrix(X)) {
    return(X)
  }

  if (is.list(X)) {
    input_names <- names(X)

    if (is.null(input_names)) {
      input_names <- as.character(seq_along(X))
    } else {
      empty <- input_names == "" | is.na(input_names)
      input_names[empty] <- as.character(which(empty))
    }

    return(Map(
      function(x, nm) RcppML::factor_input(x, nm),
      X,
      input_names
    ))
  }

  stop("X must be a matrix, sparse Matrix, or a list of such matrices.")
}

.fit_nmf_replicates <- function(
  X,
  k = 7,
  reps = 50,
  seed = NULL,
  threads = 0,
  verbose = FALSE,
  L1 = c(0.05, 0.05),
  tol = 1e-4,
  nonneg = c(TRUE, TRUE),
  ...
) {
  X <- .compose_nmf_input(X, ...)

  run_seeds <- if (!is.null(seed)) {
    withr::with_seed(seed, {
      sample.int(.Machine$integer.max, reps)
    })
  } else {
    rep(NA_integer_, reps)
  }

  fits <- vector("list", reps)

  for (r in seq_len(reps)) {
    fits[[r]] <- if (is.na(run_seeds[r])) {
      RcppML::nmf(
        X,
        k = k,
        threads = threads,
        verbose = verbose,
        L1 = L1,
        tol = tol,
        nonneg = nonneg,
        ...
      )
    } else {
      RcppML::nmf(
        X,
        k = k,
        seed = run_seeds[r],
        threads = threads,
        verbose = verbose,
        L1 = L1,
        tol = tol,
        nonneg = nonneg,
        ...
      )
    }
  }

  fits
}


.align_factors <- function(ref_fit, fit, factor_weight = 0.5) {
  ref <- .get_wh(ref_fit)
  cur <- .get_wh(fit)

  k_ref <- ncol(ref$w)
  k_cur <- ncol(cur$w)

  if (k_ref != k_cur) {
    stop(sprintf(
      "All runs must use the same k. Reference has %d, current has %d.",
      k_ref, k_cur
    ))
  }

  k <- k_ref
  sim <- matrix(0, nrow = k, ncol = k)

  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      sim_w <- .cosine_sim(ref$w[, i], cur$w[, j])
      sim_h <- .cosine_sim(ref$h[i, ], cur$h[j, ])
      sim[i, j] <- factor_weight * sim_w + (1 - factor_weight) * sim_h
    }
  }

  cost <- 1 - sim
  assignment <- as.integer(clue::solve_LSAP(cost))

  list(
    map = assignment,  # reference factor -> matched factor in this run
    similarity = sim,
    matched_similarity = sim[cbind(seq_len(k), assignment)]
  )
}




build_consensus_nmf <- function(
  X,
  k = 7,
  reps = 50,
  seed = NULL,
  threads = 0,
  verbose = FALSE,
  L1 = c(0.0, 0.0),
  factor_weight = 0.5,
  zero_tol = 1e-8,
  tol = 1e-4,
  nonneg = c(TRUE, TRUE),
  ...
) {
    
  fits <- .fit_nmf_replicates(
    X = X,
    k = k,
    reps = reps,
    seed = seed,
    threads = threads,
    verbose = verbose,
    L1 = L1,
    tol = tol,
    nonneg = nonneg,
    ...
  )

  if (is.null(fits) || length(fits) == 0) {
    stop("No fitted models were produced.")
  }

  ref_fit <- fits[[1]]
  ref_wh <- .get_wh(ref_fit)

  k <- ncol(ref_wh$w)
  n_features <- nrow(ref_wh$w)
  n_samples <- ncol(ref_wh$h)

  feature_names <- rownames(ref_wh$w)
  if (is.null(feature_names)) {
    feature_names <- paste0("feature_", seq_len(n_features))
  }

  sample_names <- colnames(ref_wh$h)
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_len(n_samples))
  }

  alignments <- vector("list", length(fits))
  alignments[[1]] <- list(
    map = seq_len(k),
    similarity = diag(k),
    matched_similarity = rep(1, k)
  )

  factor_match_similarity <- matrix(
    NA_real_,
    nrow = length(fits),
    ncol = k,
    dimnames = list(paste0("run_", seq_along(fits)), paste0("factor_", seq_len(k)))
  )
  factor_match_similarity[1, ] <- 1

  for (r in 2:length(fits)) {
    aln <- .align_factors(ref_fit, fits[[r]], factor_weight = factor_weight)
    alignments[[r]] <- aln
    factor_match_similarity[r, ] <- aln$matched_similarity
  }

  gene_votes <- matrix(
    0,
    nrow = n_features,
    ncol = k,
    dimnames = list(feature_names, paste0("factor_", seq_len(k)))
  )

  cell_votes <- matrix(
    0,
    nrow = n_samples,
    ncol = k,
    dimnames = list(sample_names, paste0("factor_", seq_len(k)))
  )

  gene_loading_sum <- matrix(
    0,
    nrow = n_features,
    ncol = k,
    dimnames = list(feature_names, paste0("factor_", seq_len(k)))
  )

  cell_loading_sum <- matrix(
    0,
    nrow = n_samples,
    ncol = k,
    dimnames = list(sample_names, paste0("factor_", seq_len(k)))
  )

  for (r in seq_along(fits)) {
    wh <- .get_wh(fits[[r]])
    map <- alignments[[r]]$map

    for (f in seq_len(k)) {
      mf <- map[f]

      gene_idx <- which(abs(wh$w[, mf]) > zero_tol)
      cell_idx <- which(abs(wh$h[mf, ]) > zero_tol)

      if (length(gene_idx)) {
        gene_votes[gene_idx, f] <- gene_votes[gene_idx, f] + 1
        gene_loading_sum[gene_idx, f] <-
          gene_loading_sum[gene_idx, f] +
          abs(wh$w[gene_idx, mf])
      }

      if (length(cell_idx)) {
        cell_votes[cell_idx, f] <- cell_votes[cell_idx, f] + 1
        cell_loading_sum[cell_idx, f] <-
          cell_loading_sum[cell_idx, f] +
          abs(wh$h[mf, cell_idx])
      }
    }
  }

  gene_support <- gene_votes / length(fits)
  cell_support <- cell_votes / length(fits)

  gene_loading_mean <- gene_loading_sum / pmax(gene_votes, 1)
  cell_loading_mean <- cell_loading_sum / pmax(cell_votes, 1)

  list(
    models = fits,
    alignments = alignments,
    factor_match_similarity = factor_match_similarity,
    gene_votes = gene_votes,
    cell_votes = cell_votes,
    gene_support = gene_support,
    cell_support = cell_support,
    gene_loading_mean = gene_loading_mean,
    cell_loading_mean = cell_loading_mean,
    feature_names = feature_names,
    sample_names = sample_names,
    k = k,
    n_runs = length(fits),
    factor_weight = factor_weight,
    zero_tol = zero_tol
  )
  }