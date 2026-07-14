bicluster_edgeR <- function(
    bic,
    counts,
    col_data,
    group_by,
    design,
    contrast,
    gene_slot = "default",
    cell_slot = "test",
    test = c("QLF", "LRT"),
    verbose = FALSE
){

    test <- match.arg(test)

    ## ------------------------------------------------------------
    ## single bicluster
    ## ------------------------------------------------------------

    is_single <-
        is.list(bic) && ("cells" %in% names(bic)) 

    if (is_single) {

        return(
            .bicluster_edgeR(
                bic = bic,
                counts = counts,
                col_data = col_data,
                group_by = group_by,
                design = design,
                contrast = contrast,
                gene_slot = gene_slot,
                cell_slot = cell_slot,
                test = test,
                verbose = verbose
            )
        )

    }

    ## ------------------------------------------------------------
    ## list of biclusters
    ## ------------------------------------------------------------

    res <- lapply(
        bic,
        .bicluster_edgeR,
        counts = counts,
        col_data = col_data,
        group_by = group_by,
        design = design,
        contrast = contrast,
        gene_slot = gene_slot,
        cell_slot = cell_slot,
        test = test,
        verbose = verbose
    )

    names(res) <- names(bic)

    res
}


.bicluster_edgeR <- function(
    bic,
    counts,
    col_data,
    group_by,
    design,
    contrast,
    test = c("QLF", "LRT"),
    cell_slot = c("test", "cells")
){
    test <- match.arg(test)

    cell_slot <- match.arg(cell_slot)

    ## ------------------------------------------------------------
    ## subset cells
    ## ------------------------------------------------------------

    cell_names <- intersect(
        names(bic[[cell_slot]]),
        colnames(counts)
    )

    counts <- counts[, cell_names, drop = FALSE]
    meta <- as.data.frame(col_data[cell_names, , drop = FALSE])

    ## ------------------------------------------------------------
    ## sanity checks
    ## ------------------------------------------------------------

    if(!all(group_by %in% colnames(meta)))
        stop("Not all group_by variables found in col_data.")

    ## ------------------------------------------------------------
    ## pseudobulk IDs
    ## ------------------------------------------------------------

    pb <- interaction(
        meta[, group_by, drop = FALSE],
        drop = TRUE
    )

    ## ------------------------------------------------------------
    ## pseudobulk counts
    ## ------------------------------------------------------------

    pb_levels <- levels(pb)

    pb_counts <- vapply(
        pb_levels,
        function(g){

            Matrix::rowSums(
                counts[, pb == g, drop = FALSE]
            )

        },
        numeric(nrow(counts))
    )

    rownames(pb_counts) <- rownames(counts)
    colnames(pb_counts) <- pb_levels

    ## ------------------------------------------------------------
    ## pseudobulk metadata
    ## ------------------------------------------------------------

    pb_meta <- unique(
        meta[, group_by, drop = FALSE]
    )

    pb_meta$.pb <- interaction(
        pb_meta[, group_by, drop = FALSE],
        drop = TRUE
    )

    rownames(pb_meta) <- pb_meta$.pb

    pb_meta <- pb_meta[pb_levels, group_by, drop = FALSE]

    pb_meta[] <- lapply(pb_meta, factor)

    ## ------------------------------------------------------------
    ## remove empty pseudobulks
    ## ------------------------------------------------------------

    keep <- Matrix::colSums(pb_counts) > 0

    if(any(!keep)){

        warning(
            sprintf(
                "Dropping %d empty pseudobulks.",
                sum(!keep)
            ),
            call. = FALSE
        )

        pb_counts <- pb_counts[, keep, drop = FALSE]
        pb_meta <- pb_meta[keep, , drop = FALSE]
    }

    ## ------------------------------------------------------------
    ## edgeR
    ## ------------------------------------------------------------
    y <- edgeR::DGEList(pb_counts)
    y <- edgeR::calcNormFactors(y)
    design_matrix <- model.matrix(
        design,
        data = pb_meta
    )

    print(dim(pb_counts))
    print(dim(design_matrix))

    y <- edgeR::estimateDisp(y, design_matrix)
    contr <- limma::makeContrasts(
        contrasts = contrast,
        levels = design_matrix
    )
    if (test == "QLF") {
        fit <- edgeR::glmQLFit(
            y,
            design_matrix
        )
        out <- edgeR::glmQLFTest(
            fit,
            contrast = contr
        )
    } else {
        fit <- edgeR::glmFit(
            y,
            design_matrix
        )
        out <- edgeR::glmLRT(
            fit,
            contrast = contr
        )
    }
    edgeR::topTags(out, n = Inf)$table
}