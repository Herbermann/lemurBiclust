.get_factors <- function(
  result,
  gene_support_threshold = 0.9,
  cell_support_threshold = 0.9
) {

  support_object <- result$consensus

  has_test <- !is.null(support_object$cell_test_support)

  gene_support <- support_object$gene_support
  cell_support <- support_object$cell_support
  if (has_test){
    test_support <- support_object$cell_test_support
    test_names <- support_object$test_names
  }


  feature_info <- .parse_feature_names(support_object$feature_names)
  sample_names <- support_object$sample_names
  
  k <- support_object$metadata$k

  factor_match_similarity <-
    support_object$metadata$factor_match_similarity

  factors <- vector("list", k)


  for (f in seq_len(k)) {

    gene_idx <- which(gene_support[, f] >= gene_support_threshold)
    cell_idx <- which(cell_support[, f] >= cell_support_threshold)

    if (has_test) {
      test_idx <- which(test_support[, f] >= cell_support_threshold)
    }

    selected_features <- feature_info[gene_idx, , drop = FALSE]

    cell_loadings <- setNames(
      support_object$cell_loading_mean[cell_idx, f],
      sample_names[cell_idx]
    )

    cell_test_loadings <- setNames(
      support_object$cell_test_loading_mean[test_idx, f],
      test_names[test_idx]
    )

    loading_by_view <- lapply(
      split(
        seq_len(nrow(selected_features)),
        ifelse(
          is.na(selected_features$view),
          "default",
          selected_features$view
        )
      ),
      function(idx) {
        setNames(
          support_object$gene_loading_mean[gene_idx[idx], f] *
            ifelse(selected_features$direction[idx] == "down", -1, 1),
          selected_features$feature[idx]
        )
      }
    )

    fac <- .init_factor()

    cells <- .init_group()
    cells$loadings <- cell_loadings
    cells$elements <- names(cell_loadings)

    test <- .init_group()
    test$loadings <- cell_test_loadings
    test$elements <- names(cell_test_loadings)

    views <- lapply(loading_by_view, function(x) {
      view <- .init_group()
      view$loadings <- x
      view$elements <- names(x)
      view
    })

    fac$views <- views
    fac$cells <- cells
    fac$test <- test
    
    fac$metadata <- list(
      id = sprintf("F%03d", f),
      match_similarity = factor_match_similarity[, f]
    )

    factors[[f]] <- fac

  }

  result$factors <- factors
  result$metadata$gene_support_threshold <- gene_support_threshold
  result$metadata$cell_support_threshold <- cell_support_threshold

  result
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

    loading_by_view <- lapply(
      split(seq_len(nrow(selected_features)),
            ifelse(is.na(selected_features$view),
                  "default",
                  selected_features$view)),
      function(idx) {
        setNames(
          support_object$gene_loading_mean[gene_idx[idx], f],
          selected_features$feature[idx]
        )
      }
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

      loading_by_view = loading_by_view,

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