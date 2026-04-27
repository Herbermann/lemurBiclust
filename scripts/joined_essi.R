library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(Matrix)
library(lemur)
library(tidyverse)
library(Polychrome)
library(scico)
library(patchwork)
library("progressr")

# Choose folders for data and results etc.
folder_data <- "/Users/herbermann/code/data/coculture/"
file_name <- "coculture_effect_in_AML_5000"


####################################################################
# Load the sce object
###################################################################

sce_all <- readRDS(paste0(folder_data, "sce_object_", file_name, ".RDS"))

###############################################################
# Run lemur
###############################################################

set.seed(42)

n_latents <- 20

fit <- lemur(sce_all, ~coculture, n_embedding=n_latents, test_fraction= 0.5)
# Storing 50% of the data (38858 cells) as test data.
# Regress out global effects using linear method.
# Find base point for differential embedding
# Fit differential embedding model
# Initial error: 1.3e+07
# ---Fit Grassmann linear model
# Final error: 8.31e+06

#fit <- align_harmony(fit)
fit <- align_by_grouping(fit, grouping = colData(fit)[,"clusters_renamed"])

fit <- test_de(fit, contrast = cond(coculture = "Coculture") - cond(coculture = "BM_only"))
nei <- find_de_neighborhoods(fit, group_by = vars(individual, coculture))
#Find optimal neighborhood using zscore.
#Validate neighborhoods using test data
#Form pseudobulk (summing counts)
#Calculate size factors for each gene
#Fit glmGamPoi model on pseudobulk data
#Fit diff-in-diff effect

pval_threshold = 0.05

nei <- nei[order(nei[,"adj_pval"]),]
nei_sig <- nei[nei[,"adj_pval"] < pval_threshold,]
dim(nei_sig)

###############################################################
# Visualize cell types etc. on the embedding
###############################################################

umap <- uwot::umap(t(fit$embedding))

what_to_visualize <- "clusters_renamed"
nbr_colors <- length(table(fit$colData[,what_to_visualize]))
#set.seed(122)
custom_colors <- createPalette(nbr_colors, seedcolors = palette36.colors(10))
#custom_colors <- sample(custom_colors)
names(custom_colors) <- levels(fit$colData[,what_to_visualize])

as_tibble(fit$colData) %>%
  mutate(umap = umap) %>%
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = .data[[what_to_visualize]]), size = 1) +
  scale_color_manual(values=custom_colors) +
  facet_wrap(vars(coculture)) +
  labs(title = "UMAP of latent space from LEMUR")


what_to_visualize <- "individual"
nbr_colors <- length(table(fit$colData[,what_to_visualize]))
custom_colors <- createPalette(nbr_colors, seedcolors = palette36.colors(10))
names(custom_colors) <- levels(fit$colData[,what_to_visualize])

as_tibble(fit$colData) %>%
  mutate(umap = umap) %>%
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = .data[[what_to_visualize]]), size = 1) +
  scale_color_manual(values=custom_colors) +
  facet_wrap(vars(coculture)) +
  labs(title = "UMAP of latent space from LEMUR")


###############################################################
# Visualize DE
###############################################################

#rownames(nei_sig) <- nei_sig[,"name"]

sel_gene <- "GATA2"
#nei_sig[sel_gene,]

#neighborhood_coordinates <- nei_sig %>%
#  dplyr::filter(name == sel_gene) %>%
#  unnest(c(neighborhood)) %>%
#  dplyr::rename(cell_id = neighborhood) %>%
#  left_join(tibble(cell_id = rownames(umap), umap), by = "cell_id") %>%
#  dplyr::select(name, cell_id, umap)

p1 <- tibble(umap = umap) %>%
  mutate(de = assay(fit, "DE")[sel_gene,]) %>%
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = de)) +
  scale_color_scico(palette="vik", midpoint=0) +
#  geom_density2d(data = neighborhood_coordinates, breaks = 0.2, 
#                 contour_var = "ndensity", color = "black") +
  labs(title = "Differential expression with neighborhood boundary")


as_tibble(fit$colData) %>%
  mutate(umap = umap) %>%
  mutate(exprs = assay(fit, "logcounts")[sel_gene,]) %>%
  ggplot(aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = exprs), size=0.5) +
  scale_color_scico(palette="vik", midpoint=0) +
  facet_wrap(vars(coculture)) +
  labs(title = sel_gene)



############################################
############################################
############################################
############################################
############################################
library(Matrix)
library("RcppML")
package.version("RcppML")


de_matrix <- as.matrix(assay(fit$training_data, "DE"))

de_pos <- pmax(de_matrix, 0.)   # genes x cells, upregulation only
de_neg <- pmax(-de_matrix, 0.)  # genes x cells, downregulation only

quant_pos <- quantile(de_pos, c(0.5, 0.75, 0.9, 0.95))
quant_neg <- quantile(de_neg, c(0.5, 0.75, 0.9, 0.95))

de_pos <- tanh(2/quant_pos[4] * de_pos)
de_neg <- tanh(2/quant_neg[4] * de_neg)

de_pos <- Matrix(de_pos, sparse=TRUE)
de_neg <- Matrix(de_neg, sparse=TRUE)




pos_orig <- rownames(de_pos)
neg_orig <- rownames(de_neg)

rownames(de_pos) <- paste0("up_", make.unique(pos_orig))
rownames(de_neg) <- paste0("down_", make.unique(neg_orig))

de <- rbind(de_pos, de_neg)



##################################

genespace_stability <- function(fit1, fit2) {
  # weight W columns by their factor importance or not?? Makes less important factors less important for overlap measure.. Think more about it.
  W1 <- fit1@w #%*% diag(fit1@d)
  W2 <- fit2@w #%*% diag(fit2@d)
  
  # remove dead factors
  W1 <- W1[, colSums(W1) > 0, drop = FALSE]
  W2 <- W2[, colSums(W2) > 0, drop = FALSE]
  
  # now QR on the weighted W
  Q1 <- qr.Q(qr(W1))
  Q2 <- qr.Q(qr(W2))
  
  decomp <- svd(t(Q1) %*% Q2)
  theta2 <- acos(decomp@d)^2
  sqrt(theta2)/(length(theta2)* 0.5 * 3.1415)
}



find_factors_old <- function(
    matrix,
    ks,
    n_restarts = 2,
    test_fraction = 0.,
    L1 = c(0., 0.),
    tol = 1e-4,
    show_plot = FALSE,
    stability_thrsh = 0.9 
){
  
  p <- progressr::progressor(steps = length(ks) * n_restarts)

  results <- lapply(ks, function(k) {
    
    # single loop over restarts, compute both error and stability from same fits
    if (test_fraction == 0.){
      fits <- lapply(seq_len(n_restarts), function(i) {
      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))
      nmf(matrix, k = k, seed = i, L1 = L1, tol = tol)
      })
    }
    else {
      fits <- lapply(seq_len(n_restarts), function(i) {
      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))
      nmf(matrix, k = k, seed = i, L1 = L1, tol = tol, test_fraction = test_fraction, mask="zeros")
      })
    }
    
    # --- reconstruction error ---
    loss_slot <- if (test_fraction == 0) "loss" else "test_loss"
    errors <- sapply(fits, function(f) f@misc[[loss_slot]])
    best_idx <- which.min(errors)
    
    # --- stability ---
    ws <- lapply(fits, function(m) {
      w <- m@w
      w[, colSums(w) > 0, drop = FALSE]
    })

    sims <- combn(length(fits), 2, function(idx) {
      if (ncol(fits[[idx[1]]]@w) == 0 || ncol(fits[[idx[2]]]@w) == 0) return(NA_real_)
      d <- genespace_stability(fits[[idx[1]]], fits[[idx[2]]])
      mean(d)
    }, simplify = TRUE)

    out <- list(
      k        = k,
      error    = errors[best_idx],
      stability = max(sims, na.rm = TRUE)
      #fit      = fits[[best_idx]]
    )
    out
  })


  if (show_plot == TRUE){

      error_df <- data.frame(
      k     = ks,
      error = sapply(results, `[[`, "error")
    )
    stab_df <- data.frame(
      k         = ks,
      stability = sapply(results, `[[`, "stability")
    )

    p_error <- ggplot(error_df, aes(x = k, y = error)) +
      geom_line() + geom_point() +
      labs(x = "k", y = "reconstruction error") +
      theme_classic()

    p_stab <- ggplot(stab_df, aes(x = k, y = stability)) +
      geom_line() + geom_point() +
      labs(x = "k", y = "mean factor stability") +
      theme_classic()

    print(p_error + p_stab)

  }

  results
}


find_factors_medium <- function(
    matrix,
    ks,
    n_restarts      = 2,
    test_fraction   = 0.,
    L1              = c(0., 0.),
    tol             = 1e-4,
    show_plot       = FALSE,
    stability_thrsh = 0.9
) {
  p <- progressr::progressor(steps = length(ks) * n_restarts)

  # --- build speckle mask once, shared across all k and restarts -----------
  # Matches RcppML's test_fraction + mask="zeros" logic:
  #   eligible pool = non-zero entries; sample test_fraction of them.
  use_cv <- test_fraction > 0
  if (use_cv) {
    nz_idx <- which(matrix != 0, arr.ind = TRUE)   # non-zero (row, col) pairs
    n_test <- round(nrow(nz_idx) * test_fraction)
    set.seed(42)                                    # reproducible mask; expose as arg if desired
    test_rows <- sample.int(nrow(nz_idx), n_test)
    test_coords <- nz_idx[test_rows, , drop = FALSE]  # matrix of (row, col) pairs
    
    # original values at test positions (for MSE computation)
    test_vals <- matrix[test_coords]
    
    # training matrix: zero out the held-out entries
    train_matrix <- matrix
    train_matrix[test_coords] <- 0
    
    test_mask <- Matrix::sparseMatrix(
      i = test_coords[, 1],
      j = test_coords[, 2],
      x = 1,
      dims = dim(matrix)
    )
  } else {
    train_matrix <- matrix
  }


  results <- lapply(ks, function(k) {

    fits <- lapply(seq_len(n_restarts), function(i) {
      nmf(train_matrix, k = k, seed = i, L1 = L1, tol = tol)
      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))
    })

    # --- reconstruction / test error ----------------------------------------
    if (use_cv) {
      # For each fit, compute MSE on the held-out entries only.
      # Prediction: W %*% diag(d) %*% H, evaluated at test_coords.
        errors <- sapply(fits, function(f) {
          evaluate(f, data = matrix, mask = test_mask, missing_only = TRUE)
        })
    } else {
      errors <- sapply(fits, function(f) f@misc[["loss"]])
    }

    best_idx <- which.min(errors)

    # --- stability -----------------------------------------------------------
    if (length(fits) < 2) {# We need to catch complications for sims when n_restarts = 1. Just make it raise NA which can be handled gracefully downstream
      sims <- NA_real_
    } else {
      sims <- combn(length(fits), 2, function(idx) {
        tryCatch({
          if (ncol(fits[[idx[1]]]@w) == 0 || ncol(fits[[idx[2]]]@w) == 0) return(NA_real_)
          mean(genespace_stability(fits[[idx[1]]], fits[[idx[2]]]))
        }, error = function(e) NA_real_)
      }, simplify = TRUE)
    }

    list(
      k         = k,
      error     = errors[best_idx],
      stability = mean(sims, na.rm = TRUE)
    )
  })

  if (show_plot) {
    error_df <- data.frame(k = ks, error     = sapply(results, `[[`, "error"))
    stab_df  <- data.frame(k = ks, stability = sapply(results, `[[`, "stability"))

    p_error <- ggplot(error_df, aes(x = k, y = error)) +
      geom_line() + geom_point() +
      labs(x = "k", y = if (use_cv) "test MSE" else "reconstruction error") +
      theme_classic()

    p_stab <- ggplot(stab_df, aes(x = k, y = stability)) +
      geom_line() + geom_point() +
      labs(x = "k", y = "mean factor stability") +
      theme_classic()

    print(p_error + p_stab)
  }

  results
}



find_factors <- function(
    matrix,
    ks,
    n_restarts      = 2,
    test_fraction   = 0.,
    L1              = c(0., 0.),
    tol             = 1e-4,
    show_plot       = FALSE,
    stability_thrsh = 0.9
) {
  p <- progressr::progressor(steps = length(ks) * n_restarts)

  use_cv <- test_fraction > 0
  if (use_cv) {
    # eligible pool computed once; mask resampled per restart
    nz_idx <- which(matrix != 0, arr.ind = TRUE)
    n_test <- round(nrow(nz_idx) * test_fraction)
  }

  results <- lapply(ks, function(k) {

    fits_and_errors <- lapply(seq_len(n_restarts), function(i) {

      if (use_cv) {
        # fresh mask per restart
        test_coords  <- nz_idx[sample.int(nrow(nz_idx), n_test), , drop = FALSE]
        test_mask    <- Matrix::sparseMatrix(
          i    = test_coords[, 1],
          j    = test_coords[, 2],
          x    = 1,
          dims = dim(matrix)
        )
        train_matrix <- matrix
        train_matrix[test_coords] <- 0
        if (inherits(matrix, "sparseMatrix"))
          train_matrix <- methods::as(train_matrix, "dgCMatrix")
      } else {
        train_matrix <- matrix
      }

      fit <- nmf(train_matrix, k = k, seed = i, L1 = L1, tol = tol)

      err <- if (use_cv) {
        evaluate(fit, data = matrix, mask = test_mask, missing_only = TRUE)
      } else {
        fit@misc[["loss"]]
      }

      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))

      list(fit = fit, error = err)
    })

    fits   <- lapply(fits_and_errors, `[[`, "fit")
    errors <- sapply(fits_and_errors, `[[`, "error")
    best_idx <- which.min(errors)

    # --- stability -----------------------------------------------------------
    if (length(fits) < 2) {
      sims <- NA_real_
    } else {
      sims <- combn(length(fits), 2, function(idx) {
        tryCatch({
          if (ncol(fits[[idx[1]]]@w) == 0 || ncol(fits[[idx[2]]]@w) == 0) return(NA_real_)
          mean(genespace_stability(fits[[idx[1]]], fits[[idx[2]]]))
        }, error = function(e) NA_real_)
      }, simplify = TRUE)
    }

    list(
      k         = k,
      error     = mean(errors),   # average over restarts, not just best
      stability = max(sims, na.rm = TRUE)
    )
  })

  if (show_plot) {
    error_df <- data.frame(k = ks, error     = sapply(results, `[[`, "error"))
    stab_df  <- data.frame(k = ks, stability = sapply(results, `[[`, "stability"))

    p_error <- ggplot(error_df, aes(x = k, y = error)) +
      geom_line() + geom_point() +
      labs(x = "k", y = if (use_cv) "test MSE" else "reconstruction error") +
      theme_classic()

    p_stab <- ggplot(stab_df, aes(x = k, y = stability)) +
      geom_line() + geom_point() +
      labs(x = "k", y = "mean factor stability") +
      theme_classic()

    print(p_error + p_stab)
  }

  results
}




with_progress({
  test <- find_factors(de, c(6, 12, 18, 20, 22), n_restarts=2, test_fraction = 0.1 , show_plot = TRUE, L1 = c(0.0 , 0.0) , tol = 1e-2)
})

# k = 15!

biclusters <- nmf(de, k = 15, seed = 1)

W <- biclusters@w
H <- biclusters@h

n_genes <- length(pos_orig)
W_pos <- W[1:n_genes, ]
W_neg  <- W[(n_genes+1):nrow(de), ]

rownames(W_pos) <- pos_orig
rownames(W_neg) <- neg_orig

hist(H[1,])
hist(W_pos[,1])
hist(W_neg[,1])



generate_null_dist <- function(matrix, k, L1 = c(0, 0), n_perm = 20,
                               thr_h = 0.95, thr_w_pos = 0.95, thr_w_neg = 0.95) {

  n_genes <- nrow(matrix) / 2

  nulls <- lapply(seq_len(n_perm), function(i) {
    mat_perm <- t(apply(matrix, 1, sample))
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-3)
    list(h = fit_perm@h, w = fit_perm@w)
  })

  null_h <- lapply(nulls, `[[`, "h")
  null_w <- lapply(nulls, `[[`, "w")

  null_w_pos <- lapply(null_w, function(w) w[1:n_genes, , drop = FALSE])
  null_w_neg <- lapply(null_w, function(w) w[(n_genes + 1):nrow(w), , drop = FALSE])

  null_threshold_h <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_h, function(h) h[i, ])), thr_h)
  })

  null_threshold_w_pos <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w_pos, function(w) w[, i])), thr_w_pos)
  })

  null_threshold_w_neg <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w_neg, function(w) w[, i])), thr_w_neg)
  })

  list(h = null_threshold_h, w_pos = null_threshold_w_pos, w_neg = null_threshold_w_neg)
}

null_dist <- generate_null_dist(de, k = 15, n_perm = 5, thr_h = c(0.9, 0.95, 0.99, 0.995, 0.999), thr_w_pos = c(0.9, 0.95, 0.99, 0.995, 0.999), thr_w_neg = c(0.9, 0.95, 0.99, 0.995, 0.999))

# decide which CL to use!
null_h <- null_dist$h[3,]
null_w_pos <- null_dist$w_pos[3,]
null_w_neg <- null_dist$w_neg[3,]


n_bic <- nrow(H)

bic_cells <- lapply(1:n_bic, function(i){
  H[i,] > null_h[i]
})

bic_up <- lapply(1:n_bic, function(i){
  W_pos[,i] > null_w_pos[i]
})

bic_down <- lapply(1:n_bic, function(i){
  W_neg[,i] > null_w_neg[i]
})




#################
# Plottings.
#umap <- uwot::umap(t(fit$training_data$embedding))

umap <-uwot::umap(t(assay(fit$training_data,"DE")))

##############
k <- 3
l <- 1

de_expression <- de_matrix

genes_up <- sort(W_pos[,k], decreasing = TRUE)

genes_down <- sort(W_neg[,k], decreasing = TRUE)

sel_gene <- names(genes_up[l])
#sel_gene <- "GATA2"



p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit$training_data, "DE")[sel_gene,]) |>
  #arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()



p2 <- tibble(umap = umap) |>
  mutate(member_cells = bic_cells[[k]]) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = member_cells), size = 0.5) +
    coord_fixed()

p1 | p2



################
jaccard_with_logic <- function(A, B){
  sum(A & B) / sum(A | B)
}

compute_jaccard_matrices <- function(W, H, null_h, null_w) {
  
  k <- nrow(H)
  
  # binary memberships for all factors
  cell_members <- lapply(seq_len(k), function(i) H[i, ] > null_h[i])
  gene_members <- lapply(seq_len(k), function(i) W[, i] > null_w[i])
  
  # k x k jaccard matrices
  jaccard_cells <- outer(seq_len(k), seq_len(k), Vectorize(function(i, j) {
    jaccard_with_logic(cell_members[[i]], cell_members[[j]])
  }))
  
  jaccard_genes <- outer(seq_len(k), seq_len(k), Vectorize(function(i, j) {
    jaccard_with_logic(gene_members[[i]], gene_members[[j]])
  }))
  
  rownames(jaccard_cells) <- colnames(jaccard_cells) <- paste0("k", seq_len(k))
  rownames(jaccard_genes) <- colnames(jaccard_genes) <- paste0("k", seq_len(k))
  
  list(cells = jaccard_cells, genes = jaccard_genes)
}

jac_up <- compute_jaccard_matrices(W_pos, H, null_h, null_w_pos)
jac_up$cells
jac_up$genes


jac_down <- compute_jaccard_matrices(W_neg, H, null_h, null_w_neg)
jac_down$cells
jac_down$genes

pheatmap(jac_up$cells)
pheatmap(jac_up$genes)
pheatmap(jac_down$genes)

######
#Investigate connection of biclusters to known annotations
celltype_map <- setNames(sce_all$clusters_renamed, colnames(sce_all))
patient_map <- setNames(sce_all$individual, colnames(sce_all))

bic_celltypes <- lapply(1:n_bic, function(i){
  celltype_map[which(bic_cells[[i]] == TRUE)]
})

plot_df <- bind_rows(lapply(seq_along(bic_celltypes), function(i) {
  ct <- bic_celltypes[[i]]
  
  df <- as.data.frame(table(ct))
  colnames(df) <- c("celltype", "n")
  
  df$bic <- paste0("Bic ", i)
  df
}))

#Investigate connection of biclusters to known annotations
# --- align all objects to the same cells ---
common_cells <- Reduce(intersect, c(
  list(names(celltype_map), names(patient_map)),
  lapply(bic_cells, names)
))

celltype_map2 <- celltype_map[common_cells]
patient_map2   <- patient_map[common_cells]
bic_cells2     <- lapply(bic_cells, function(x) x[common_cells])

# --- cell type composition within each bicluster ---
bic_celltypes <- lapply(seq_along(bic_cells2), function(i) {
  celltype_map2[bic_cells2[[i]]]
})

plot_df <- bind_rows(lapply(seq_along(bic_celltypes), function(i) {
  ct <- bic_celltypes[[i]]

  df <- as.data.frame(table(ct))
  colnames(df) <- c("celltype", "n")
  df$bic <- paste0("Bic ", i)
  df
})) %>%
  filter(!is.na(celltype)) %>%
  group_by(bic) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(plot_df, aes(x = bic, y = prop, fill = celltype)) +
  geom_col(width = 0.8) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Fraction of cells in bicluster",
    fill = "Cell type",
    title = "Cell type composition across biclusters"
  ) +
  theme_minimal(base_size = 13)

# --- optional: patient composition within each bicluster ---
bic_patient <- lapply(seq_along(bic_cells2), function(i) {
  patient_map2[bic_cells2[[i]]]
})

plot_df_patient <- bind_rows(lapply(seq_along(bic_patient), function(i) {
  df <- as.data.frame(table(patient = bic_patient[[i]]))
  colnames(df) <- c("patient", "n")
  df$bic <- paste0("Bic ", i)
  df
})) %>%
  filter(!is.na(patient)) %>%
  group_by(bic) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(plot_df_patient, aes(x = bic, y = prop, fill = patient)) +
  geom_col(width = 0.8) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Fraction of cells in bicluster",
    fill = "Patient",
    title = "Patient composition across biclusters"
  ) +
  theme_minimal(base_size = 13)

############

score_bicluster <- function(de_matrix, membership, standardize = TRUE, eps = 1e-8) {
  stopifnot(ncol(de_matrix) == length(membership))

  in_idx  <- which(membership)
  out_idx <- which(!membership)

  if (length(in_idx) == 0) stop("membership contains no selected cells.")
  if (length(out_idx) == 0) stop("membership selects all cells; no background remains.")

  x_in  <- de_matrix[, in_idx, drop = FALSE]
  x_out <- de_matrix[, out_idx, drop = FALSE]

  score <- rowMeans(x_in, na.rm = TRUE) - rowMeans(x_out, na.rm = TRUE)

  if (standardize) {
    s_in  <- matrixStats::rowSds(as.matrix(x_in),  na.rm = TRUE)
    s_out <- matrixStats::rowSds(as.matrix(x_out), na.rm = TRUE)
    pooled <- sqrt((s_in^2 + s_out^2) / 2)
    score <- score / (pooled + eps)
  }

  names(score) <- rownames(de_matrix)
  sort(score, decreasing = TRUE)
}


k_ <- 1
bic_up <- W_pos[,k_] > null_w_pos[k_]
bic_down <- W_neg[,k_] > null_w_neg[k_]
bic_cells <- H[k_,] > null_h[k_]
bic_score <- score_bicluster(de_matrix, bic_cells)


library(fgsea)
library(msigdbr)


scores <- score_bicluster(de_matrix, bic_cells)

msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

res <- fgsea(pathways = pathways, stats = scores)
res <- res[order(res$padj), ]
head(res, 10)


scores_list <- lapply(1:ncol(W),
  function(i) {
    score_bicluster(de_matrix, (H[i,] > null_h[i]))
})

res_list <- lapply(1:ncol(W),
  function(i) {
    fgsea(pathways = pathways, stats = scores_list[[i]])
})

bic_names <- paste0("BIC", seq_along(res_list))
names(res_list) <- bic_names

all_pathways <- unique(unlist(lapply(res_list, function(x) x$pathway)))

nes_mat <- matrix(
  NA_real_,
  nrow = length(res_list),
  ncol = length(all_pathways),
  dimnames = list(bic_names, all_pathways)
)

for (i in seq_along(res_list)) {
  res <- res_list[[i]]
  nes_mat[i, res$pathway] <- res$NES
}

library(pheatmap)

pheatmap(
  nes_mat,
  scale = "none", 
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
d <- dist(nes_mat)
hc <- hclust(d)

plot(hc)

cor_mat <- cor(t(nes_mat))
pheatmap(cor_mat)


table <- do.call(
  rbind,
  lapply(seq_along(res_list), function(i) {
    g <- res_list[[i]]

    if (is.null(g) || nrow(g) == 0) return(NULL)
    
    g_sig <- g[g$padj < 0.05, ]
    if (nrow(g_sig) == 0) return(NULL)  # <- important
    
    data.frame(
      bicluster = i,   # <- rename for clarity
      pathway = g_sig$pathway,
      NES = g_sig$NES,
      padj = g_sig$padj,
      pval = g_sig$pval,
      ES = g_sig$ES,
      size = g_sig$size,
      stringsAsFactors = FALSE
    )
  })
)
table

pca_nes <- prcomp(nes_mat,center = TRUE, scale. = TRUE)
plot(pca_nes$x[,1], pca_nes$x[,2])

pca_nes$rotation[sort(abs(pca_nes$rotation[,4]), decreasing = TRUE),4]
summary(pca_nes)


gata2 <- de_matrix["GATA2",]
summary(gata2)

gata2_bic <- sapply(seq_along(bic_cells), function(i) {
  mean(gata2[bic_cells[[i]]], na.rm = TRUE)
})

names(gata2_bic) <- paste0("BIC", seq_along(bic_cells))
gata2_bic

cor(gata2_bic, pca_nes$x)
pc1 <- pca_nes$x[,1]
pc1
plot_df$bic <- gsub(" ", "", plot_df$bic)  # "Bic 1" → "Bic1"
plot_df$bic <- toupper(plot_df$bic)        # "Bic1" → "BIC1"
plot_df$pc1 <- pc1[plot_df$bic]
celltype_pc1 <- aggregate(
  prop * pc1 ~ celltype,
  data = plot_df,
  FUN = sum
)

colnames(celltype_pc1)[2] <- "PC1_score"
celltype_pc1
barplot(
  celltype_pc1$PC1_score,
  names.arg = celltype_pc1$celltype,
  las = 2,
  main = "Cell types along PC1",
  ylab = "PC1 (activation axis)"
)
plot(pca_nes$x[,1], gata2_bic, pch=19)
