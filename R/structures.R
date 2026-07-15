
# --- group ---
# representation of of elements of a bicluster
# used for cells, test cells, different gene views etc.
.init_group <- function() {
  structure(
    list(
      elements = character(0),
      loadings = numeric(0),
      metadata = list()
    ),
    class = "grp"
  )
}


# --- factor ---
# representation of both factors and biclusters
.init_factor <- function() {
  structure(
    list(
      views = list(),
      cells = .init_group(),
      test  = NULL,
      metadata = list()
    ),
    class = "biclust_factor"
  )
}


# --- BiclustResult ---
# container class for Biclustering result
# further functions act mostly on this container class

.init_BiclustResult <- function() {
  structure(
    list(
      consensus       = list(),
      factors         = list(),
      biclusters      = list(),
      summary         = list(),
      diagnostics     = list(),
      analyses        = list(),
      metadata        = list()
    ),
    class = "BiclustResult"
  )
}


#' @export
print.BiclustResult <- function(x, ...) {

    cat("\n")
    cat("<BiclustResult>\n\n")

    cat(sprintf(
        " Consensus factors    : %d\n",
        length(x$factors)
    ))

    cat(sprintf(
        " Candidate biclusters : %d\n",
        length(x$biclusters)
    ))

    if (length(x$analyses) == 0) {

        cat(" Analyses             : none\n")

    } else {

        cat(sprintf(
            " Analyses             : %s\n",
            paste(names(x$analyses), collapse = ", ")
        ))

    }

    cat("\nAccessors:\n")
    cat("  biclusters()\n")
    cat("  factors()\n")
    cat("  analyses()\n")
    cat("  diagnostics()\n")
    cat("  summary()\n")

    invisible(x)
}

#' @export
summary.BiclustResult <- function(object, ...) {

    cat("\n")
    cat("<BiclustResult summary>\n")
    cat("-------------------------\n\n")

    cat(sprintf(
        "Consensus factors      : %d\n",
        length(object$factors)
    ))

    cat(sprintf(
        "Candidate biclusters   : %d\n",
        length(object$biclusters)
    ))

    if (!is.null(object$metadata$k)) {

        cat(sprintf(
            "Requested factors (k)  : %d\n",
            object$metadata$k
        ))

    }

    if (!is.null(object$metadata$gene_support_threshold)) {

        cat(sprintf(
            "Gene support threshold : %.2f\n",
            object$metadata$gene_support_threshold
        ))

    }

    if (!is.null(object$metadata$cell_support_threshold)) {

        cat(sprintf(
            "Cell support threshold : %.2f\n",
            object$metadata$cell_support_threshold
        ))

    }

    ## views
    if (length(object$biclusters) > 0) {

        bc <- biclusters(object, 1)

        views <- setdiff(
            names(bc),
            c("cells", "test")
        )

        cat(sprintf(
            "Views                  : %s\n",
            paste(views, collapse = ", ")
        ))

    }

    ## analyses
    cat("\nAnalyses\n")

    if (length(object$analyses) == 0) {

        cat("  none\n")

    } else {

        for (nm in names(object$analyses)) {

            cat(sprintf(
                "  ✓ %s\n",
                nm
            ))

        }

    }

    invisible(object)
}

# --- consensus ---
# container for storing a consensus NMF run
# contains all details about individual runs,
# as well as matching between different runsand further
.init_consensus <- function() {
    structure(
        list(
            models = list(),
            alignments = NULL,

            gene_support = NULL,
            cell_support = NULL,

            gene_loading_mean = NULL,
            cell_loading_mean = NULL,

            sample_names = NULL,
            feature_names = NULL,

            metadata = list()
        ),
        class = "biclust_consensus"
    )
}