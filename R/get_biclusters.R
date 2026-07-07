.parse_feature_names <- function(x) {

  parts <- strsplit(x, "|", fixed = TRUE)

  info <- lapply(seq_along(parts), function(i) {

    p <- parts[[i]]

    if (length(p) >= 3 && p[2] %in% c("up", "down")) {

      data.frame(
        id = x[i],
        view = p[1],
        direction = p[2],
        feature = paste(p[-c(1, 2)], collapse = "|"),
        stringsAsFactors = FALSE
      )

    } else if (length(p) >= 2 && p[1] %in% c("up", "down")) {

      data.frame(
        id = x[i],
        view = NA_character_,
        direction = p[1],
        feature = paste(p[-1], collapse = "|"),
        stringsAsFactors = FALSE
      )

    } else {

      data.frame(
        id = x[i],
        view = NA_character_,
        direction = NA_character_,
        feature = paste(p, collapse = "|"),
        stringsAsFactors = FALSE
      )

    }
  })

  do.call(rbind, info)
}


get_biclusters <- function(
  support_object,
  gene_support_threshold = 0.9,
  cell_support_threshold = 0.9
) {

  has_test <- !is.null(support_object$cell_test_support)

  gene_support <- support_object$gene_support
  cell_support <- support_object$cell_support
  if (has_test){
    test_support <- support_object$cell_test_support
  }

  feature_info <- .parse_feature_names(support_object$feature_names)
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

  if (has_test) {
    bicluster_summary$n_test <- integer(k)
  }
  
  for (f in seq_len(k)) {

    gene_idx <- which(gene_support[, f] >= gene_support_threshold)
    cell_idx <- which(cell_support[, f] >= cell_support_threshold)

    if (has_test) {
      test_idx <- which(test_support[, f] >= cell_support_threshold)
    }

    selected_features <- feature_info[gene_idx, , drop = FALSE]

    gene_loadings <- setNames(
      support_object$gene_loading_mean[gene_idx, f],
      selected_features$feature
    )

    cell_loadings <- setNames(
      support_object$cell_loading_mean[cell_idx, f],
      sample_names[cell_idx]
    )

    bc <- list(

      factor = f,
      feature_info = selected_features,

      gene_names = unique(selected_features$feature),

      gene_up = unique(
        selected_features$feature[
          selected_features$direction == "up"
        ]
      ),

      gene_down = unique(
        selected_features$feature[
          selected_features$direction == "down"
        ]
      ),

      genes_by_view = split(
        selected_features$feature,
        ifelse(
          is.na(selected_features$view),
          "default",
          selected_features$view
        )
      ),

      cell_names = sample_names[cell_idx],

      gene_loading = gene_loadings,
      cell_loading = cell_loadings
    )

    if (has_test) {

      test_names <- support_object$test_names

      bc$test_names <- test_names[test_idx]

      bc$test_loading <- setNames(
        support_object$cell_test_loading_mean[test_idx, f],
        test_names[test_idx]
      )

      bicluster_summary$n_test[f] <- length(test_idx)
    }

    biclusters[[f]] <- bc

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
  }

  list(
    biclusters = biclusters,
    summary = bicluster_summary,
    gene_support_threshold = gene_support_threshold,
    cell_support_threshold = cell_support_threshold
  )
}