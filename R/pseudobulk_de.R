.pseudobulk_bic <- function(lemur_test, count_assay, cells, genes, design, group_vars = NULL) {
  count_assay <- SummarizedExperiment::assay(lemur_test, count_assay)

  cells <- base::intersect(as.character(cells), colnames(count_assay))
  genes <- base::intersect(as.character(genes), rownames(count_assay))

  if (length(cells) == 0 || length(genes) == 0) {
    stop("No genes or cells for bic.")
  }

  meta <- as.data.frame(lemur_test$colData)
  meta$cell_id <- rownames(meta)

  md <- meta[match(cells, meta$cell_id), , drop = FALSE]
  if (anyNA(md$cell_id)) {
    stop("Some cells were not found in the metadata.")
  }

  if (is.null(group_vars)) {
    group_vars <- all.vars(delete.response(terms(design)))
  }

  if (!all(group_vars %in% names(md))) {
    stop(
      "Some group_vars are not present in metadata: ",
      paste(setdiff(group_vars, names(md)), collapse = ", ")
    )
  }

  md[group_vars] <- lapply(md[group_vars], function(x) {
    if (!is.factor(x)) factor(x) else x
  })

  md$pb_id <- interaction(md[, group_vars, drop = FALSE],
                          drop = TRUE, sep = "_", lex.order = TRUE)

  mat <- count_assay[genes, md$cell_id, drop = FALSE]
  pb <- t(rowsum(t(mat), group = md$pb_id, reorder = FALSE))

  sample_table <- md[!duplicated(md$pb_id), c("pb_id", group_vars), drop = FALSE]
  rownames(sample_table) <- sample_table$pb_id
  sample_table <- sample_table[colnames(pb), , drop = FALSE]

  list(
    pb = pb,
    sample_table = sample_table,
    design = design
  )
}



# We wanna make it understand lemur based condition contrasts!
.as_edger_contrast <- function(contrast) {
  if (inherits(contrast, "contrast_relation")) {
    contrast <- evaluate_contrast_tree(
      contrast,
      contrast,
      function(x, y) x
    )
  }

  if (inherits(contrast, "model_vec")) {
    contrast <- unclass(contrast)
  }

  contrast
}



.run_edger <- function(pb_counts, sample_table, design, contrast) {
  sample_table <- as.data.frame(sample_table)
  sample_table <- sample_table[colnames(pb_counts), , drop = FALSE]
  sample_table[] <- lapply(sample_table, function(x) if (is.character(x)) factor(x) else x)

  design_pb <- model.matrix(design, data = sample_table)

  y <- edgeR::DGEList(counts = pb_counts)
  keep <- edgeR::filterByExpr(y, design = design_pb)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  y <- edgeR::estimateDisp(y, design_pb)
  fit <- edgeR::glmQLFit(y, design_pb)

  if (is.null(contrast)){
    cntrst <- evaluate_contrast_tree(
      lemur_fit$contrast,
      lemur_fit$contrast,
      function(x, y) x - y
    )
  }

  contrast <- .as_edger_contrast(contrast)

  stopifnot(any(contrast != 0))

  res <- edgeR::glmQLFTest(fit, contrast = contrast)

  edgeR::topTags(res, n = Inf)$table
}