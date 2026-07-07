.pareto_core <- function(loadings,
                        fallback_explained = 0.10,
                        min_keep = 20) {
  x <- loadings[is.finite(loadings) & loadings > 0]
  x <- sort(x, decreasing = TRUE)

  n <- length(x)
  if (n == 0) {
    return(list(
      core_cells = character(0),
      core_idx = integer(0),
      method = "none",
      core_fraction = NA_real_,
      loading_explained = NA_real_,
      elbow_idx = NA_integer_,
      elbow_explained = NA_real_
    ))
  }

  frac_cells <- seq_len(n) / n
  frac_loading <- cumsum(x) / sum(x)

  gain <- frac_loading - frac_cells
  elbow_idx <- which.max(gain)
  elbow_explained <- frac_loading[elbow_idx]

  if (elbow_explained >= fallback_explained) {
    k <- elbow_idx
    method <- "elbow"
  } else {
    k <- which(frac_loading >= fallback_explained)[1]
    if (is.na(k)) k <- n
    method <- "coverage"
  }

  k <- max(k, min_keep)
  k <- min(k, n)

  list(
    selected = names(x)[seq_len(k)],
    core_idx = seq_len(k),
    method = method,
    core_fraction = unname(k / n),
    loading_explained = unname(frac_loading[k]),
    elbow_idx = unname(elbow_idx),
    elbow_explained = unname(elbow_explained)
  )

}


.extract_core_bicluster <- function(
  bc,
  core_method = c("pareto", "kmeans")
) {

  core_method <- match.arg(core_method)

  if (core_method == "kmeans") {
    stop("Not implemented")
  }

  out_gene <- .pareto_core(bc$gene_loading)
  out_cell <- .pareto_core(bc$cell_loading)

  core <- list(
    genes = out_gene,
    cells = out_cell
  )

  if (!is.null(bc$test_loading)) {
    core$test <- .pareto_core(bc$test_loading)
  }

  bc$core <- core

  bc
}


extract_core <- function(
  result_object,
  core_method = c("pareto", "kmeans")
) {

  core_method <- match.arg(core_method)

  biclusters <- lapply(
    result_object$biclusters,
    .extract_core_bicluster,
    core_method = core_method
  )

  has_test <- !is.null(biclusters[[1]]$core$test)

  k <- length(biclusters)

  summary <- data.frame(
    bicluster = seq_len(k),
    n_genes = integer(k),
    n_cells = integer(k),
    stringsAsFactors = FALSE
  )

  if (has_test)
    summary$n_test <- integer(k)

  for (i in seq_len(k)) {

    summary$n_genes[i] <-
      length(biclusters[[i]]$core$genes$selected)

    summary$n_cells[i] <-
      length(biclusters[[i]]$core$cells$selected)

    if (has_test) {
      summary$n_test[i] <-
        length(biclusters[[i]]$core$test$selected)
    }
  }

  result_object$biclusters <- biclusters
  result_object$core_summary <- summary

  result_object
}