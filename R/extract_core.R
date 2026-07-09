.pareto_core <- function(
  in_group,
  fallback_explained = 0.10,
  min_keep = 20
) {

  signed_loadings <- in_group$loadings

  keep <- is.finite(signed_loadings) & signed_loadings != 0

  x <- signed_loadings[keep]
  ord <- order(abs(x), decreasing = TRUE)

  x <- x[ord]

  out_group <- .init_group()

  n <- length(x)

  if (n == 0) {

    out_group$metadata <- list(
      method = "none",
      core_fraction = NA_real_,
      loading_explained = NA_real_,
      elbow_idx = NA_integer_,
      elbow_explained = NA_real_
    )

    return(out_group)
  }

  frac_cells <- seq_len(n) / n
  frac_loading <- cumsum(abs(x)) / sum(abs(x))

  gain <- frac_loading - frac_cells
  elbow_idx <- which.max(gain)
  elbow_explained <- frac_loading[elbow_idx]

  if (elbow_explained >= fallback_explained) {
    k <- elbow_idx
    method <- "elbow"
  } else {
    k <- which(frac_loading >= fallback_explained)[1]
    if (is.na(k)) {
      k <- n
    }
    method <- "coverage"
  }

  k <- max(k, min_keep)
  k <- min(k, n)

  out_group$elements <- names(x)[seq_len(k)]
  out_group$loadings <- x[seq_len(k)]

  out_group$metadata <- list(
    method = method,
    core_fraction = unname(k / n),
    loading_explained = unname(frac_loading[k]),
    elbow_idx = unname(elbow_idx),
    elbow_explained = unname(elbow_explained)
  )

  out_group
}



.extract_bicluster <- function(
  fac,
  core_method = c("pareto", "kmeans")
) {

  core_method <- match.arg(core_method)

  if (core_method == "kmeans") {
    stop("Not implemented")
  }

  bc <- .init_factor()

  bc$cells <- .pareto_core(fac$cells)

  bc$views <- lapply(
    fac$views,
    .pareto_core
  )

  if (!is.null(fac$test)) {
    bc$test <- .pareto_core(fac$test)
  }

  bc$metadata <- fac$metadata

  bc
}


.extract_core <- function(
  result,
  core_method = c("pareto", "kmeans")
) {

  core_method <- match.arg(core_method)

  factors <- result$factors

  biclusters <- lapply(
    factors,
    .extract_bicluster,
    core_method = core_method
  )

  result$biclusters <- biclusters

  result
}





