

get_biclusters <- function(
  support_object,
  assay = NULL,
  gene_support_threshold = 0.50,
  cell_support_threshold = 0.50,
  null_validate = NULL
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
      gene_names = sub("^(up|down)_", "", feature_names[gene_idx]),
      gene_up = sub("^up_", "", feature_names[gene_idx][grepl("^up_", feature_names[gene_idx])]),
      gene_down = sub("^down_", "", feature_names[gene_idx][grepl("^down_", feature_names[gene_idx])]),
      cell_names = sample_names[cell_idx]
    )

    if (!is.null(assay)) {
      gs <- biclusters[[f]]$gene_names
      cs <- biclusters[[f]]$cell_names

      if (length(gs) < 2 || length(cs) < 2) {
        coh <- NA_real_
      } else {
        sub <- assay[gs, cs, drop = FALSE]
        coh <- irlba::irlba(sub, nv = 1, nu = 0)$d^2 / norm(sub, type = "F")^2
      }

      if (is.numeric(null_validate)) {
        global_coherence <- irlba::irlba(assay, nv = 1, nu = 0)$d^2 / norm(assay, type = "F")^2
        reps <- null_validate
        null_coherences <- numeric(reps)

        for (j in seq_len(reps)) {
          if (length(gs) < 2 || length(cs) < 2) {
            null_coherences[j] <- NA_real_
          } else {
            row_idx <- sample.int(nrow(assay), length(gs), replace = FALSE)
            col_idx <- sample.int(ncol(assay), length(cs), replace = FALSE)
            Xnull <- assay[row_idx, col_idx, drop = FALSE]
            null_coherences[j] <- irlba::irlba(Xnull, nv = 1, nu = 0)$d^2 / norm(Xnull, type = "F")^2
          }
        }

        sd <- sqrt(mean((null_coherences - global_coherence)^2, na.rm = TRUE))
        p_emp <- (sum(null_coherences >= coh, na.rm = TRUE) + 1) / (sum(!is.na(null_coherences)) + 1)
      } else {
        sd <- NULL
        p_emp <- NULL
        global_coherence <- NULL
      }

    } else {
      coh <- NULL
      sd <- NULL
      p_emp <- NULL
      global_coherence <- NULL
    }

    bicluster_summary$n_genes[f] <- length(gene_idx)
    bicluster_summary$n_cells[f] <- length(cell_idx)
    bicluster_summary$mean_gene_support[f] <- mean(gene_support[, f], na.rm = TRUE)
    bicluster_summary$mean_cell_support[f] <- mean(cell_support[, f], na.rm = TRUE)
    bicluster_summary$max_gene_support[f] <- max(gene_support[, f], na.rm = TRUE)
    bicluster_summary$max_cell_support[f] <- max(cell_support[, f], na.rm = TRUE)
    bicluster_summary$mean_match_similarity[f] <- mean(factor_match_similarity[, f], na.rm = TRUE)
    bicluster_summary$coherence[f] <- coh
    bicluster_summary$null_coherence[f] <- global_coherence
    bicluster_summary$null_deviation[f] <- sd
    bicluster_summary$p_emp[f] <- p_emp
  }

  list(
    biclusters = biclusters,
    summary = bicluster_summary,
    gene_support_threshold = gene_support_threshold,
    cell_support_threshold = cell_support_threshold
  )
}

