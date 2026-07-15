#' @export
biclusters<- function(result, bic = NULL) {

  if (!inherits(result, "BiclustResult"))
    stop(
        "'result' must be a BiclustResult.",
        call. = FALSE
    )

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
  structure(
    out,
    class = "BiclusterList"
  )
}





#' @export
analyses <- function(result) {

    if (!inherits(result, "BiclustResult"))
        stop("'result' must be a BiclustResult.",
             call. = FALSE)

    result$analyses
}


#' @export
consensus <- function(result) {

    if (!inherits(result, "BiclustResult"))
        stop("'result' must be a BiclustResult.",
             call. = FALSE)

    result$consensus
}


#' @export
factors <- function(result) {

    if (!inherits(result, "BiclustResult"))
        stop("'result' must be a BiclustResult.",
             call. = FALSE)

    result$factors
}


#' @export
diagnostics <- function(result) {

    if (!inherits(result, "BiclustResult"))
        stop("'result' must be a BiclustResult.",
             call. = FALSE)

    result$diagnostics
}






#' Cell composition of biclusters
#'
#' Summarize the composition of biclusters with respect to one or more
#' metadata variables.
#'
#' @param bic A bicluster or list of biclusters.
#' @param x SingleCellExperiment or lemur object.
#' @param label Character vector giving one or more colData columns.
#' @param cell_slot Which cells to use.
#'
#' @export
bicluster_composition <- function(
    bic,
    x,
    label,
    cell_slot = c("test", "cells", "both")
){

    cell_slot <- match.arg(cell_slot)

    meta <- as.data.frame(SummarizedExperiment::colData(x))

    is_single <-
        is.list(bic) &&
        ("cells" %in% names(bic))

    if(is_single){

        out <- .bicluster_composition(
            bic,
            meta,
            label = label,
            cell_slot = cell_slot
        )

        out$bicluster <-
            if(!is.null(bic$metadata$id))
                bic$metadata$id
            else
                "B001"

        return(out)

    }

    out <- lapply(
        seq_along(bic),
        function(i){

            res <- .bicluster_composition(
                bic[[i]],
                meta,
                label = label,
                cell_slot = cell_slot
            )

            res$bicluster <-
                if(!is.null(names(bic)))
                    names(bic)[i]
                else
                    sprintf("B%03d", i)

            res

        }
    )

    dplyr::bind_rows(out)

}


.bicluster_composition <- function(
    bic,
    meta,
    label,
    cell_slot
){

    slots <-

        switch(
            cell_slot,
            test = "test",
            cells = "cells",
            both = c("cells", "test")
        )

    out <- lapply(
        slots,
        function(slot){

            cells <- intersect(
                names(bic[[slot]]),
                rownames(meta)
            )

            df <- meta[cells, label, drop = FALSE]

            res <- lapply(
                label,
                function(var){

                    n_total <- table(meta[[var]])

                    n_bic <- table(df[[var]])

                    lev <- union(
                        names(n_total),
                        names(n_bic)
                    )

                    n_total <- n_total[lev]
                    n_bic <- n_bic[lev]

                    n_total[is.na(n_total)] <- 0
                    n_bic[is.na(n_bic)] <- 0

                    data.frame(
                        variable = var,
                        label = lev,
                        n = as.integer(n_bic),
                        bic_size = length(cells),
                        label_size = as.integer(n_total),
                        fraction_bic =
                            as.integer(n_bic) /
                            length(cells),
                        fraction_label =
                            as.integer(n_bic) /
                            as.integer(n_total),
                        stringsAsFactors = FALSE
                    )

                }
            )

            res <- dplyr::bind_rows(res)

            res$cell_slot <- slot

            res

        }
    )

    dplyr::bind_rows(out)

}