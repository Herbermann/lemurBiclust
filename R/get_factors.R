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
      cell_test_loadings <- setNames(
        support_object$cell_test_loading_mean[test_idx, f],
        test_names[test_idx]
      )}

    selected_features <- feature_info[gene_idx, , drop = FALSE]

    cell_loadings <- setNames(
      support_object$cell_loading_mean[cell_idx, f],
      sample_names[cell_idx]
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
        
        sign <- ifelse(
          selected_features$direction[idx] == "down", -1, 1)
        sign[is.na(sign)] <- 1
          
        x <- setNames(
          support_object$gene_loading_mean[gene_idx[idx], f] * sign,
          selected_features$feature[idx]
        )
        ## Collapse duplicated gene flavours by keeping the
        ## loading with the largest absolute value.
        if (anyDuplicated(names(x))) {
          x <- tapply(
            x,
            names(x),
            function(w) w[which.max(abs(w))]
          )
          x <- unlist(x, use.names = TRUE)
        }
        x
      }
    )

    fac <- .init_factor()

    cells <- .init_group()
    cells$loadings <- cell_loadings
    cells$elements <- names(cell_loadings)

    fac$cells <- cells

    if (has_test) {

      cell_test_loadings <- setNames(
        support_object$cell_test_loading_mean[test_idx, f],
        test_names[test_idx]
      )

      test <- .init_group()
      test$loadings <- cell_test_loadings
      test$elements <- names(cell_test_loadings)

      fac$test <- test
    }

    views <- lapply(loading_by_view, function(x) {
      view <- .init_group()
      view$loadings <- x
      view$elements <- names(x)
      view
    })

    fac$views <- views
    
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


