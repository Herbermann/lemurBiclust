#' @export
bicluster_edgeR <- function(x, ...) {
    UseMethod("bicluster_edgeR")
}


#' @export
bicluster_edgeR.BiclusterList <- function(
    x,
    lemur_fit,
    use_assay = "counts",
    group_by,
    design,
    contrast,
    gene_slot = "default",
    cell_slot = "test",
    test = c("QLF", "LRT"),
    verbose = FALSE
){
    bic <- x
  
    sce <- switch(
        cell_slot,
        test = lemur_fit$test_data,
        stop("Only 'test' is currently supported.")
    )

    counts   <- SummarizedExperiment::assay(sce, use_assay)
    col_data <- SummarizedExperiment::colData(sce)

    ## ------------------------------------------------------------
    ## Prepare formula for LEMUR contrast parsing
    ## ------------------------------------------------------------

    mf <- stats::model.frame(design, data = as.data.frame(col_data))

    mm <- stats::model.matrix(design, data = mf)

    attr(design, "xlevels") <- stats::.getXlevels(stats::terms(design), mf)
    attr(design, "vars_xlevels") <- lapply(mf, function(x) if (is.factor(x)) levels(x) else NULL)
    attr(design, "contrasts") <-  attr(mm, "contrasts")

    contrast <- .parse_contrast(
        contrast = {{ contrast }},
        formula = design
    )

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


#' @export
bicluster_edgeR.BiclustResult <- function(
    x,
    lemur_fit,
    use_assay = "counts",
    group_by,
    design,
    contrast,
    gene_slot = "default",
    cell_slot = "test",
    test = c("QLF", "LRT"),
    verbose = FALSE
){
    result <- x
  
    sce <- switch(
        cell_slot,
        test = lemur_fit$test_data,
        stop("Only 'test' is currently supported.")
    )

    counts   <- SummarizedExperiment::assay(sce, use_assay)
    col_data <- SummarizedExperiment::colData(sce)

    ## ------------------------------------------------------------
    ## Prepare formula for LEMUR contrast parsing
    ## ------------------------------------------------------------

    mf <- stats::model.frame(design, data = as.data.frame(col_data))
    mm <- stats::model.matrix(design, data = mf)

    attr(design, "xlevels") <- stats::.getXlevels(stats::terms(design), mf)
    attr(design, "vars_xlevels") <- lapply(mf, function(x) if (is.factor(x)) levels(x) else NULL)
    attr(design, "contrasts") <-  attr(mm, "contrasts")

    contrast <- .parse_contrast(
        contrast = {{ contrast }},
        formula = design
    )


    bic <- biclusters(result)

    result$analyses$edgeR <- bicluster_edgeR(
        x = bic,
        lemur_fit = lemur_fit,
        use_assay = use_assay,
        group_by = group_by,
        design = design,
        contrast = contrast,
        gene_slot = gene_slot,
        cell_slot = cell_slot,
        test = test,
        verbose = verbose
    )

    result
}


.bicluster_edgeR <- function(
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

    if (!cell_slot %in% names(bic))
      stop(sprintf("cell_slot '%s' not found.", cell_slot))

    if (!gene_slot %in% names(bic))
        stop(sprintf("gene_slot '%s' not found.", gene_slot))
      
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

    y <- edgeR::estimateDisp(y, design_matrix)

    if (test == "QLF") {
        fit <- edgeR::glmQLFit(y, design_matrix)
        out <- edgeR::glmQLFTest(
            fit,
            contrast = contrast
        )
    } else {
        fit <- edgeR::glmFit(y, design_matrix)
        out <- edgeR::glmLRT(
            fit,
            contrast = contrast
        )
    }

    tab <- edgeR::topTags(out, n = Inf)$table
    tab$gene <- rownames(tab)

    tab
}