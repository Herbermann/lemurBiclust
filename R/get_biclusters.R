
get_biclusters <- function(
  support_object,
  assay = NULL,
  gene_support_threshold = 0.9,
  cell_support_threshold = 0.9
) {


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

    gene_loadings <- setNames(
      support_object$gene_loading_mean[gene_idx, f],
      sub("^(up|down)_", "", feature_names[gene_idx])
    )

    cell_loadings <- setNames(
      support_object$cell_loading_mean[cell_idx, f],
      sample_names[cell_idx]
    )

    biclusters[[f]] <- list(
      factor = f,

      gene_names = sub(
        "^(up|down)_",
        "",
        feature_names[gene_idx]
      ),

      gene_up = sub(
        "^up_",
        "",
        feature_names[gene_idx][grepl("^up_", feature_names[gene_idx])]
      ),

      gene_down = sub(
        "^down_",
        "",
        feature_names[gene_idx][grepl("^down_", feature_names[gene_idx])]
      ),

      cell_names = sample_names[cell_idx],

      gene_loading = gene_loadings,
      cell_loading = cell_loadings
    )

    if (!is.null(assay)) {

      gs <- biclusters[[f]]$gene_names
      cs <- biclusters[[f]]$cell_names

      if (length(gs) < 2 || length(cs) < 2) {

        coh <- NA_real_

      } else {

        sub <- assay[gs, cs, drop = FALSE]

        coh <- irlba::irlba(
          sub,
          nv = 1,
          nu = 0
        )$d^2 / norm(sub, type = "F")^2
      }

    } else {

      coh <- NULL
    }

    bicluster_summary$n_genes[f] <- length(gene_idx)
    bicluster_summary$n_cells[f] <- length(cell_idx)

    bicluster_summary$mean_gene_support[f] <-
      mean(gene_support[, f], na.rm = TRUE)

    bicluster_summary$mean_cell_support[f] <-
      mean(cell_support[, f], na.rm = TRUE)

    bicluster_summary$max_gene_support[f] <-
      max(gene_support[, f], na.rm = TRUE)

    bicluster_summary$max_cell_support[f] <-
      max(cell_support[, f], na.rm = TRUE)

    bicluster_summary$mean_match_similarity[f] <-
      mean(factor_match_similarity[, f], na.rm = TRUE)

    bicluster_summary$coherence[f] <- coh
  }

  list(
    biclusters = biclusters,
    summary = bicluster_summary,
    gene_support_threshold = gene_support_threshold,
    cell_support_threshold = cell_support_threshold
  )
}







get_test_biclusters <- function(
  support_object,
  assay = NULL,
  gene_support_threshold = 0.9,
  cell_support_threshold = 0.9
) {

  gene_support <- support_object$gene_support
  cell_support <- support_object$cell_support
  cell_test_support <- support_object$cell_test_support

  feature_names <- support_object$feature_names
  sample_names <- support_object$sample_names
  test_names <- support_object$test_names
  k <- support_object$k

  factor_match_similarity <- support_object$factor_match_similarity

  biclusters <- vector("list", k)

  bicluster_summary <- data.frame(
    factor = seq_len(k),
    n_genes = integer(k),
    n_cells = integer(k),
    n_test = integer(k),
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
    test_idx <- which(cell_test_support[, f] >= cell_support_threshold)

    gene_loadings <- setNames(
      support_object$gene_loading_mean[gene_idx, f],
      sub("^(up|down)_", "", feature_names[gene_idx])
    )

    cell_loadings <- setNames(
      support_object$cell_loading_mean[cell_idx, f],
      sample_names[cell_idx]
    )

    test_loadings <- setNames(
      support_object$cell_test_loading_mean[test_idx, f],
      test_names[test_idx]
    )

    biclusters[[f]] <- list(
      factor = f,

      gene_names = sub(
        "^(up|down)_",
        "",
        feature_names[gene_idx]
      ),

      gene_up = sub(
        "^up_",
        "",
        feature_names[gene_idx][grepl("^up_", feature_names[gene_idx])]
      ),

      gene_down = sub(
        "^down_",
        "",
        feature_names[gene_idx][grepl("^down_", feature_names[gene_idx])]
      ),

      cell_names = sample_names[cell_idx],
      test_names = test_names[test_idx],

      gene_loading = gene_loadings,
      cell_loading = cell_loadings,
      test_loading = test_loadings
    )

    if (!is.null(assay)) {

      gs <- biclusters[[f]]$gene_names
      cs <- biclusters[[f]]$cell_names

      if (length(gs) < 2 || length(cs) < 2) {

        coh <- NA_real_

      } else {

        sub <- assay[gs, cs, drop = FALSE]

        coh <- irlba::irlba(
          sub,
          nv = 1,
          nu = 0
        )$d^2 / norm(sub, type = "F")^2
      }

    } else {

      coh <- NULL
    }

    bicluster_summary$n_genes[f] <- length(gene_idx)
    bicluster_summary$n_cells[f] <- length(cell_idx)
    bicluster_summary$n_test[f] <- length(test_idx)


    bicluster_summary$mean_gene_support[f] <-
      mean(gene_support[, f], na.rm = TRUE)

    bicluster_summary$mean_cell_support[f] <-
      mean(cell_support[, f], na.rm = TRUE)

    bicluster_summary$max_gene_support[f] <-
      max(gene_support[, f], na.rm = TRUE)

    bicluster_summary$max_cell_support[f] <-
      max(cell_support[, f], na.rm = TRUE)

    bicluster_summary$mean_match_similarity[f] <-
      mean(factor_match_similarity[, f], na.rm = TRUE)

    bicluster_summary$coherence[f] <- coh
  }

  list(
    biclusters = biclusters,
    summary = bicluster_summary,
    gene_support_threshold = gene_support_threshold,
    cell_support_threshold = cell_support_threshold
  )
}
