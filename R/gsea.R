
#' Gene set enrichment analysis
#'
#' Perform GSEA using fgsea for one or more ranked gene lists.
#'
#' @param pathways Named list of gene sets to test enrichment of.
#' @param stats A named numeric vector or a list of named numeric vectors representing a gene ranking.
#' @param ... Passed to fgsea::fgsea().
#'
#' @export
custom_gsea <- function(
    pathways,
    stats,
    ...
){
    if(!requireNamespace("fgsea", quietly = TRUE))
        stop(
            "Please install 'fgsea'.",
            call. = FALSE
        )
    # --- single ranking ---

    if(is.numeric(stats)){
        if(is.null(names(stats)))
            stop(
                "'stats' must be a named numeric vector.",
                call. = FALSE
            )

        stats <- stats[!is.na(stats)]
        stats <- sort(stats, decreasing = TRUE)

        return(
            as.data.frame(
                fgsea::fgsea(
                    pathways = pathways,
                    stats = stats,
                    ...
                )
            )
        )
    }

    # --- list of rankings ---

    if(!is.list(stats))
        stop(
            "'stats' must be either a named numeric vector or a list thereof.",
            call. = FALSE
        )

    res <- lapply(
        seq_along(stats),
        function(i){

            x <- stats[[i]]

            if(is.null(names(x)))
                stop(
                    sprintf(
                        "Ranking %d is not named.",
                        i
                    ),
                    call. = FALSE
                )

            x <- x[!is.na(x)]
            x <- sort(x, decreasing = TRUE)

            out <- as.data.frame(
                fgsea::fgsea(
                    pathways = pathways,
                    stats = x,
                    ...
                )
            )

            out$ranking <-
                if(!is.null(names(stats)))
                    names(stats)[i]
                else
                    i

            out

        }
    )

    dplyr::bind_rows(res)

}



#' edgeR gene ranking
#' 
#' Create a gene ranking based on differential expression analysis from edgeR.
#' Computes a signed statistic.
#' 
#' @param de Output from 'bicluster_edgeR'
#' @param statistic Quantity to base the signed statistic on, default 'logFC'. Other options involve test statistics 'F' and 'LR'.
#' 
#' @export
edgeR_ranking <- function(
    de,
    statistic = c("logFC", "F", "LR")
){

    statistic <- match.arg(statistic)

    make_one <- function(x){

        if(statistic == "logFC"){

            stats <- x$logFC

        } else{

            if(!statistic %in% colnames(x))
                stop(
                    sprintf(
                        "Statistic '%s' not found.",
                        statistic
                    ),
                    call. = FALSE
                )

            stats <-
                sign(x$logFC) *
                sqrt(x[[statistic]])

        }

        if("gene" %in% colnames(x)){

            names(stats) <- x$gene

        } else{

            names(stats) <- rownames(x)

        }

        sort(stats, decreasing = TRUE)

    }

    if(is.data.frame(de))
        return(make_one(de))

    lapply(de, make_one)

}


#' Extract gene sets from biclusters
#' 
#' Convenience function to extract all gene sets from a list of biclusters
#' 
#' @param bic A single bicluster, or a list of biclusters, as returned by 'biclusters()'
#' 
#' @export
bicluster_gene_sets <- function(
    bic
) {
    .bicluster_gene_sets(bic)$out
}


.bicluster_gene_sets <- function(
    bic
){
    if ("cells" %in% names(bic)){
        is_single <- TRUE
    } else {
        is_single <- FALSE
    }

    if (is_single){
        out_ <- .single_bicluster_gene_sets(bic)
        out <- out_$out
        n_views <- out$n_views
    } else {
        out_ <- lapply(
            bic,
            function(x){.single_bicluster_gene_sets(x)}
        )
        out <- lapply(out_, function(y){y$out})
        n_views <- out_[[1]]$n_views
    }
    
    return(list("out" = out, "n_views" = n_views, "is_single" = is_single))
}



.single_bicluster_gene_sets <- function(
    bic,
    split_sign = TRUE
) {

    views <- names(bic)
    views <- views[!views %in% c("cells", "test")]

    n_views <- length(views)

    out <- setNames(
        lapply(
            views,
            function(name){

                w <- bic[[name]]

                if(!split_sign){

                    names(w)

                } else{

                    list(
                        up = names(w[w > 0]),
                        down = names(w[w < 0])
                    )

                }

            }
        ),
        views
    )

    list(
        "out" = out,
        "n_views" = n_views
    )
}




#' GSEA of bicluster gene sets
#' 
#' Perform GSEA of bicluster gene sets using fgsea.
#' 
#' @param x bicluster or list of biclusters, as returned by bicluster() or BiclustResult object
#' @param de Result of 'bicluster_edgeR()' used for creating a signed ranking of genes
#' 
#' @export
bicluster_gsea <- function(x, ...) {
    UseMethod("bicluster_gsea")
} 

#' @export
bicluster_gsea.BiclusterList <- function(
    bics,
    de
) {
    stats <- edgeR_ranking(de)
    gene_sets <- .bicluster_gene_sets(bics)
    n_views <- gene_sets$n_views
    is_single <- gene_sets$is_single
    gene_sets <- gene_sets$out

    if (is_single) {

        if (n_views == 1) {

            pathways <- gene_sets

            if (is.list(pathways[[1]])) {
                pathways <- unlist(pathways, recursive = FALSE)
            }

            temp <- custom_gsea(
                pathways = pathways,
                stats = stats
            )

        } else {

            temp <- lapply(
                names(gene_sets),
                function(view){

                    pathways <- gene_sets[[view]]

                    if (is.list(pathways[[1]])) {
                        pathways <- unlist(
                            pathways,
                            recursive = FALSE
                        )
                    }

                    out <- custom_gsea(
                        pathways = pathways,
                        stats = stats
                    )

                    out$view <- view

                    out

                }
            )

            temp <- dplyr::bind_rows(temp)

        }

    } else {

        temp <- lapply(
            seq_along(gene_sets),
            function(i){

                if (n_views == 1){

                    pathways <- gene_sets[[i]]

                    if (is.list(pathways[[1]])) {
                        pathways <- unlist(
                            pathways,
                            recursive = FALSE
                        )
                    }

                    out <- custom_gsea(
                        pathways = pathways,
                        stats = stats[[i]]
                    )

                } else {

                    out <- lapply(
                        names(gene_sets[[i]]),
                        function(view){

                            pathways <- gene_sets[[i]][[view]]

                            if (is.list(pathways[[1]])) {
                                pathways <- unlist(
                                    pathways,
                                    recursive = FALSE
                                )
                            }

                            x <- custom_gsea(
                                pathways = pathways,
                                stats = stats[[i]]
                            )

                            x$view <- view

                            x

                        }
                    )

                    out <- dplyr::bind_rows(out)

                }

                out

            }
        )
    }

    dplyr::bind_rows(temp, .id = "bicluster")
}



#' @export
bicluster_gsea.BiclustResult <- function(
    result,
    de = NULL
){

    analyses <- analyses(result)

    if (is.null(de)){
      if (is.null(analyses$edgeR)) {
          stop(
              "No edgeR analysis found. Run bicluster_edgeR() first.",
              call. = FALSE
          ) } else {
            edge <- analyses$edgeR
      }
    } else {
        edge <- de
    }
  
    result$analyses$gsea <- bicluster_gsea(
        biclusters(result),
        edge
    )

    result
}





annotate_biclusters <- function(
  bics,
  pathways
) {
  
}