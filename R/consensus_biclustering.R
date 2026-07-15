
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


.fit_nmf_replicates <- function(
  X,
  k = 7,
  reps = 50,
  seed = NULL,
  threads_RcppML = 1,
  verbose = FALSE,
  L1 = c(0.05, 0.05),
  tol = 1e-4,
  ...
) {

  run_seeds <- if (!is.null(seed)) {
    withr::with_seed(seed, {
      sample.int(.Machine$integer.max, reps)
    })
  } else {
    rep(NA_integer_, reps)
  }

  fits <- vector("list", reps)


  p <- progressr::progressor(steps = reps)


  for (r in seq_len(reps)) {
    fits[[r]] <- if (is.na(run_seeds[r])) {
      RcppML::nmf(
        X,
        k = k,
        threads = threads_RcppML,
        L1 = L1,
        tol = tol,
        ...
      )
    } else {
      RcppML::nmf(
        X,
        k = k,
        seed = run_seeds[r],
        threads = threads_RcppML,
        L1 = L1,
        tol = tol,
        ...
      )
    }

    p(message = sprintf("Restart %d/%d", r, reps))

  }

  fits
}



.fit_nmf_replicates_parallel <- function(
  X,
  k = 7,
  reps = 50,
  seed = NULL,
  workers = future::availableCores(),
  verbose = FALSE,
  L1 = c(0.05, 0.05),
  tol = 1e-4,
  ...
) {

  run_seeds <- if (!is.null(seed)) {
    withr::with_seed(seed, {
      sample.int(.Machine$integer.max, reps)
    })
  } else {
    rep(NA_integer_, reps)
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = workers)

  p <- progressr::progressor(steps = reps)

  fits <- future.apply::future_lapply(
    seq_len(reps),
    future.seed = TRUE,
    future.scheduling = 1,
    FUN = function(r) {

      fit <- if (is.na(run_seeds[r])) {
        RcppML::nmf(
          X,
          k = k,
          threads = 1,
          L1 = L1,
          tol = tol,
          ...
        )
      } else {
        RcppML::nmf(
          X,
          k = k,
          seed = run_seeds[r],
          threads = 1,
          L1 = L1,
          tol = tol,
          ...
        )
      }

      p(sprintf("Restart %d/%d", r, reps))

      fit
    }
  )

  fits
}


.fit_nmf_replicates_progress <- function(...){
  progressr::handlers("cli")
  progressr::with_progress({.fit_nmf_replicates(...)})
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


.align_replicates <- function(
  fits,
  factor_weight = 0.5
) {

  ref_fit <- fits[[1]]
  ref_wh <- .get_wh(ref_fit)

  k <- ncol(ref_wh$w)

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
    dimnames = list(
      paste0("run_", seq_along(fits)),
      paste0("factor_", seq_len(k))
    )
  )

  factor_match_similarity[1, ] <- 1

  if (length(fits) > 1) {
    for (r in 2:length(fits)) {

      aln <- .align_factors(
        ref_fit = ref_fit,
        fit = fits[[r]],
        factor_weight = factor_weight
      )

      alignments[[r]] <- aln
      factor_match_similarity[r, ] <- aln$matched_similarity
    }
  }

  list(
    alignments = alignments,
    factor_match_similarity = factor_match_similarity
  )
}


.factor_consensus <- function(
  fits,
  alignments,
  zero_tol = 1e-6
) {

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
      matched_factor <- map[f]

      gene_idx <- which(abs(wh$w[, matched_factor]) > zero_tol)
      cell_idx <- which(abs(wh$h[matched_factor, ]) > zero_tol)

      if (length(gene_idx)) {
        gene_votes[gene_idx, f] <- gene_votes[gene_idx, f] + 1
        gene_loading_sum[gene_idx, f] <-
          gene_loading_sum[gene_idx, f] +
          abs(wh$w[gene_idx, matched_factor])
      }

      if (length(cell_idx)) {
        cell_votes[cell_idx, f] <- cell_votes[cell_idx, f] + 1
        cell_loading_sum[cell_idx, f] <-
          cell_loading_sum[cell_idx, f] +
          abs(wh$h[matched_factor, cell_idx])
      }
    }
  }

  gene_support <- gene_votes / length(fits)
  cell_support <- cell_votes / length(fits)

  gene_loading_mean <- gene_loading_sum / pmax(gene_votes, 1)
  cell_loading_mean <- cell_loading_sum / pmax(cell_votes, 1)


  consensus <- .init_consensus()

  consensus$models <- fits
  consensus$alignments <- alignments

  consensus$gene_support <- gene_support
  consensus$cell_support <- cell_support

  consensus$gene_loading_mean <- gene_loading_mean
  consensus$cell_loading_mean <- cell_loading_mean

  consensus$metadata <- list(
      k = k,
      n_runs = length(fits),
      zero_tol = zero_tol
  )

  consensus$feature_names <- feature_names
  consensus$sample_names <- sample_names

  return(consensus)
}




.build_consensus_nmf <- function(
  X,
  k = 7,
  reps = 20,
  seed = NULL,
  threads = 0,
  threads_RcppML = 1,
  verbose = FALSE,
  L1 = c(0.0, 0.0),
  factor_weight = 0.5,
  zero_tol = 1e-8,
  tol = 1e-4,
  backend = c("serial", "parallel"),
  ...
) {
  backend <- match.arg(backend)

  fit_fun <- switch(
    backend,
    serial = .fit_nmf_replicates,
    parallel = .fit_nmf_replicates_parallel
  )

  fits <- progressr::with_progress({
  fit_fun(
    X = X,
    k = k,
    reps = reps,
    seed = seed,
    threads_RcppML = threads_RcppML,
    L1 = L1,
    tol = tol,
    ...
  )
  })

  if (is.null(fits) || length(fits) == 0) {
    stop("No fitted models were produced.")
  }

  aligned <- .align_replicates(
    fits,
    factor_weight = factor_weight
  )

  consensus <- .factor_consensus(
    fits,
    alignments = aligned$alignments,
    zero_tol = zero_tol
  )

  consensus$metadata$factor_match_similarity <-
    aligned$factor_match_similarity
    
  res <- .init_BiclustResult()
  
  res$consensus <- consensus
  
  res
}



#' @export
build_consensus_nmf <- function(
  X,
  k = 7,
  reps = 50,
  seed = NULL,
  threads = 0,
  threads_RcppML = 1,
  verbose = FALSE,
  L1 = c(0.0, 0.0),
  factor_weight = 0.5,
  zero_tol = 1e-6,
  tol = 1e-2,
  backend = c("serial", "parallel"),
  ...
) {

  .build_consensus_nmf(
    X,
    k = k,
    reps = reps,
    seed = seed,
    threads = threads,
    threads_RcppML = threads_RcppML,
    verbose = verbose,
    L1 = L1,
    factor_weight = factor_weight,
    zero_tol = zero_tol,
    tol = tol,
    backend = backend,
    ...
  )
}