

.pareto_core <- function(loadings,
                        fallback_explained = 0.60,
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
    core_cells = names(x)[seq_len(k)],
    core_idx = seq_len(k),
    method = method,
    core_fraction = k / n,
    loading_explained = frac_loading[k],
    elbow_idx = elbow_idx,
    elbow_explained = elbow_explained
  )

}



extract_core <- function(
  result_object,
  core_method = c("pareto", "kmeans")
) {

  core_method <- match.arg(core_method)
  core_cells_by_bic <- vector("list", length(result_object$biclusters))
  core_genes_by_bic <- vector("list", length(result_object$biclusters))

  if ("test_loading" %in% names(result_object$biclusters[[1]])){
    has_test <- TRUE
    core_cells_test_by_bic <- vector("list", length(result_object$biclusters))
    core_genes_test_by_bic <- vector("list", length(result_object$biclusters))
  } else {
    has_test <- FALSE
  }

  if (core_method == "pareto"){

    for (i in seq_along(result_object$biclusters)) {
      
      out_cells <-  .pareto_core(result_object$biclusters[[i]]$cell_loading)
      core_cells_by_bic[[i]] <- out_cells$core_cells

      out_genes <- .pareto_core(result_object$biclusters[[i]]$gene_loading)
      core_genes_by_bic[[i]] <- out_genes$core_cells

      if (has_test){
        out_cells <- .pareto_core(result_object$biclusters[[i]]$test_loading)
        core_cells_test_by_bic[[i]] <- out_cells$core_cells
      }
    }

  } else if (core_method == "kmeans"){
    stop("Not implemented")
  }

  if (has_test) {
    return(list(
      train = core_cells_by_bic,
      test = core_cells_test_by_bic,
      genes = core_genes_by_bic
    ))
  } else {
    return(list(
      train = core_cells_by_bic,
      genes = core_genes_by_bic
    ))
  }
}