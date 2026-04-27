devtools::load_all("/Users/herbermann/code/lemur_biclusters/lemur")
library("RcppML")
#######

library("lemur")
library("tidyverse")
library("SingleCellExperiment")
library("patchwork")
library("RcppML")
library("progressr")


set.seed(42)

data("glioblastoma_example_data", package = "lemur")
glioblastoma_example_data

fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition, 
             n_embedding = 15, test_fraction = 0.5)

fit <- align_harmony(fit)

umap <- uwot::umap(t(fit$embedding))

fit <- compute_contrasts(fit, contrast = cond(condition = "panobinostat") - cond(condition = "ctrl"))

sel_gene <- "ENSG00000169429" # is CXCL8

p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")[sel_gene,]) |>
  arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()


######################
# LEGACY NEIGHBORHOODS

neighborhoods <- find_de_neighborhoods(fit, group_by = vars(patient_id, condition))

as_tibble(neighborhoods) |>
  left_join(as_tibble(rowData(fit)[,1:2]), by = c("name" = "gene_id")) |>
  relocate(symbol, .before = "name") |>
  arrange(pval) |>
  head(5)


##################################
# ------ PREPARE MATRIX ---------#
##################################

library(Matrix)

de_matrix <- as.matrix(assay(fit, "DE"))

de_pos <- pmax(de_matrix, 0.2)   # genes x cells, upregulation only
de_neg <- pmax(-de_matrix, 0.2)  # genes x cells, downregulation only

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
  # weight W columns by their factor importance
  W1 <- fit1@w %*% diag(fit1@d)
  W2 <- fit2@w %*% diag(fit2@d)
  
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



find_factors <- function(
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
      stability = mean(sims, na.rm = TRUE)
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



#with_progress({
#  test <- find_factors(de, c(12, 13, 14, 15, 16), n_restarts=5, test_fraction = 0.1, show_plot = TRUE, L1 = c(0.0 , 0.0) )
#})

# k = 15!

biclusters <- nmf(de, k = 15, seed = 1:10)

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
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-4)
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

null_dist <- generate_null_dist(de, k = 15, n_perm = 10)

null_h <- null_dist$h
null_w_pos <- null_dist$w_pos
null_w_neg <- null_dist$w_neg


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


##############
k <- 12
l <- 1

de_expression <- de_matrix

genes_up <- sort(W_pos[,k], decreasing = TRUE)

genes_down <- sort(W_neg[,k], decreasing = TRUE)

sel_gene <- names(genes_up[l])
#sel_gene <- "GATA2"


p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")[sel_gene,]) |>
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

