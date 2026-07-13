#' @export
biclusters <- function(result, bic = NULL) {

  bcs <- result$biclusters
  ids <- vapply(bcs, function(x) x$metadata$id, character(1))

  if (is.null(bic)) {
    idx <- seq_along(bcs)
  } else if (is.numeric(bic)) {
    if (any(bic < 1 | bic > length(bcs) | bic != as.integer(bic))) {
      stop(
        sprintf(
          "'bic' must contain integers between 1 and %d.",
          length(bcs)
        ),
        call. = FALSE
      )
    }

    idx <- as.integer(bic)

  } else if (is.character(bic)) {

    idx <- match(bic, ids)

    if (anyNA(idx)) {
      stop(
        sprintf("Unknown bicluster ID(s): %s", paste(bic[is.na(idx)], collapse = ", ")),call. = FALSE
      )
    }
  } else {
    stop("'bic' must be NULL, an integer, a character string, or a vector thereof.", call. = FALSE)
  }

  out <- lapply(idx, function(i) {
    bc <- bcs[[i]]
    obj <- list()
    obj$cells <- bc$cells$loadings
    for (view in names(bc$views)) {
      obj[[view]] <- bc$views[[view]]$loadings
    }
    if (!is.null(bc$test)) {
      obj$test <- bc$test$loadings
    }
    obj
  })
  names(out) <- ids[idx]
  if (length(out) == 1L) {
    return(out[[1]])
  }
  out
}