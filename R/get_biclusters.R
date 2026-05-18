

get_biclusters <- function(
  support_object,
  gene_support_threshold = 0.50,
  cell_support_threshold = 0.50
) {
  required_fields <- c(
    "gene_support", "cell_support", "feature_names", "sample_names", "k", "factor_match_similarity"
  )

  missing_fields <- setdiff(required_fields, names(support_object))
  if (length(missing_fields) > 0) {
    stop(
      "support_object is missing required fields: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  gene_support <- support_object$gene_support
  cell_support <- support_object$cell_support
  feature_names <- support_object$feature_names
  sample_names <- support_object$sample_names
  k <- support_object$k
  factor_match_similarity <- support_object$factor_match_similarity

  biclusters <- vector("list", k)
  bicluster_summary <- data.frame(
    factor = seq_len(k),
    n_genes = integer(k),
    n_cells = integer(k),
    mean_gene_support = numeric(k),
    mean_cell_support = numeric(k),
    max_gene_support = numeric(k),
    max_cell_support = numeric(k),
    mean_match_similarity = numeric(k),
    stringsAsFactors = FALSE
  )

  for (f in seq_len(k)) {
    gene_idx <- which(gene_support[, f] >= gene_support_threshold)
    cell_idx <- which(cell_support[, f] >= cell_support_threshold)

    biclusters[[f]] <- list(
      factor = f,
      gene_idx = gene_idx,
      cell_idx = cell_idx,
      gene_names = sub("^(up|down)_", "", feature_names[gene_idx]),
      cell_names = sample_names[cell_idx],
      gene_support = gene_support[gene_idx, f, drop = FALSE],
      cell_support = cell_support[cell_idx, f, drop = FALSE]
    )

    bicluster_summary$n_genes[f] <- length(gene_idx)
    bicluster_summary$n_cells[f] <- length(cell_idx)
    bicluster_summary$mean_gene_support[f] <- mean(gene_support[, f], na.rm = TRUE)
    bicluster_summary$mean_cell_support[f] <- mean(cell_support[, f], na.rm = TRUE)
    bicluster_summary$max_gene_support[f] <- max(gene_support[, f], na.rm = TRUE)
    bicluster_summary$max_cell_support[f] <- max(cell_support[, f], na.rm = TRUE)
    bicluster_summary$mean_match_similarity[f] <- mean(factor_match_similarity[, f], na.rm = TRUE)
  }

  list(
    biclusters = biclusters,
    summary = bicluster_summary,
    gene_support_threshold = gene_support_threshold,
    cell_support_threshold = cell_support_threshold
  )
}

