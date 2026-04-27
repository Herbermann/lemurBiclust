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

fit <- lemur(sce_all, ~coculture, n_embedding=n_latents, test_fraction=0.6)
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





##############


genespace_stability <- function(W1, W2) {
  # principal angles between column spaces of W1 and W2
  # high similarity = same subspace even if individual factors differ
  Q1 <- qr.Q(qr(W1))
  Q2 <- qr.Q(qr(W2))
  decomp <- svd(t(Q1) %*% Q2)
  #mean(decomp@d)
  theta2 <- acos(decomp@d)^2
  sqrt(theta2)/(length(theta2)*0.5 * 3.1415)
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
        nmf(matrix, k = k, seed = i, L1 = L1, tol = tol, test_fraction = test_fraction)
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

    sims <- combn(length(ws), 2, function(idx) {
      W1 <- ws[[idx[1]]]
      W2 <- ws[[idx[2]]]
      if (ncol(W1) == 0 || ncol(W2) == 0) return(NA_real_)
      #sim_mat <- cosine(W1, W2)
      sim_mat <- genespace_stability(W1,W2)
      #mean(apply(sim_mat, 1, max))
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




sweep_l1 <- function(matrix,
                     k,
                     l1_values,
                     n_restarts = 2,
                     test_fraction = 0.1,
                     tol = 1e-4) {
  
  p <- progressr::progressor(steps = length(l1_values) * n_restarts)
  
  results <- lapply(l1_values, function(l1) {
    fits <- lapply(seq_len(n_restarts), function(i) {
      p(message = sprintf("L1=%.3f, restart %d/%d", l1, i, n_restarts))
      nmf(matrix, k = k, seed = i, L1 = c(l1, 0), tol = tol)
    })
    errors <- sapply(fits, function(f) f@misc$loss)
    best <- fits[[which.min(errors)]]@misc$loss
    
    list(
      l1        = l1,
      test_loss = best
    )
  })
  
  df <- data.frame(
    l1        = sapply(results, `[[`, "l1"),
    test_loss = sapply(results, `[[`, "test_loss")
  )
  
  p <- ggplot(df, aes(x = l1, y = test_loss)) +
    geom_line() + geom_point() +
    labs(x = "L1 (W)", y = "test loss", title = sprintf("L1 sweep at k=%d", k)) +
    theme_classic()
  
  list(plot = p, data = df)
}






with_progress({
  test <- find_factors(de_pos, c(12, 15, 17, 18), n_restarts=2, show_plot = TRUE, L1 = c(0.05 , 0.0) )
})



with_progress({
  test_L <- sweep_l1(de_pos, k = 12, l1_values = c(0.00001, 0.001))
})
test_L$plot


#biclusters <- nmf(de_pos, 12, L1 = c(0., 0.), seed = 1:5)
#biclusters_lasso <- nmf(de_pos, 18, L1 = c(0.05, 0.))
#biclusters_elast <- nmf(de_pos, 12, L1 = c(0.05, 0.), L2=c(0.001, 0.))



#bl1 <- nmf(de_pos, k=12, L1 = c(0.05, 0.))
#bl2 <- nmf(de_pos, k=12, L1 = c(0.1, 0.))
#bl3 <- nmf(de_pos, k=12, L1 = c(0.075, 0.))


mean_overlap <- function(fit, top_n = 50) {
  W <- fit@w
  # for each factor, get indices of top n genes by loading
  top_genes <- apply(W, 2, function(w) order(w, decreasing = TRUE)[seq_len(top_n)])
  
  # pairwise Jaccard similarity of top gene sets
  pairs <- combn(ncol(W), 2, function(idx) {
    g1 <- top_genes[, idx[1]]
    g2 <- top_genes[, idx[2]]
    length(intersect(g1, g2)) / length(union(g1, g2))
  }, simplify = TRUE)
  
  mean(pairs)
}

###################


biclusters <- nmf(de_pos, k = 12, seed = 1:5)
W <- biclusters@w
H <- biclusters@h

rownames(W) <- rownames(de_matrix)
colnames(H) <- colnames(de_matrix)

hist(W[,12])
hist(H[,12])


###############
generate_null_dist <- function(matrix, k, L1 = c(0, 0), n_perm = 20, thr_h=0.95, thr_w=0.95) {
  
  nulls <- lapply(seq_len(n_perm), function(i) {
    mat_perm <- t(apply(matrix, 1, sample))
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-2)
    list(h = fit_perm@h, w = fit_perm@w)
  })

  null_h <- lapply(nulls, `[[`, "h")
  null_w <- lapply(nulls, `[[`, "w")

  null_threshold_h <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_h, function(h) h[i, ])), thr_h)
  })

  null_threshold_w <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w, function(w) w[, i])), thr_w)
  })

  list(h = null_threshold_h, w = null_threshold_w)
}

null_pos_ <- generate_null_dist(de_pos, k = 12, n_perm = 10, thr_h = c(0.95,0.99), thr_w = c(0.95,0.99))

null_pos   <- null_pos_$h[2,]
null_pos_W <- null_pos_$w[2,]


#################
# Plottings.
umap <- uwot::umap(t(fit$training_data$embedding))

umap <-uwot::umap(t(assay(fit$training_data,"DE")))


k <- 5
l <- 1

de_expression <- de_matrix

genes <- sort(W[,k], decreasing = TRUE)
genes[l]
sel_gene <- names(genes[l])
#sel_gene <- "GATA2"


p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit$training_data, "DE")[sel_gene,]) |>
  #arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()


null_thresh_h <- null_pos[k]
null_thresh_w <- null_pos_W[k]

member_cells <- H[k,]>null_pos[k]
member_genes <- W[,k]>null_pos_W[k]


p2 <- tibble(umap = umap) |>
  mutate(member_cells = member_cells) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = member_cells), size = 0.5) +
    coord_fixed()

p1 | p2


out_nei <- de_expression[sel_gene, which(member_cells == FALSE)]
in_nei <- de_expression[sel_gene, which(member_cells == TRUE)]

t.test(in_nei, out_nei)



num_bic <- length(W[1,])

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

jac <- compute_jaccard_matrices(W, H, null_pos, null_pos_W)
jac$cells
jac$genes

extract_genes <- function(W, k){ 
  names(sort(W[,k], decreasing = TRUE)>null_pos_W[k])
}

extract_genes(W, k)[1:40]

############

score_bicluster <- function(de_matrix, membership, standardize = TRUE, eps = 1e-8) {
  stopifnot(ncol(de_matrix) == length(membership))

  cells_in <- which(membership)

  if (length(cells_in) == 0) {
    stop("membership contains no selected cells.")
  }
  if (length(cells_in) == ncol(de_matrix)) {
    stop("membership selects all cells; no background remains.")
  }

  x <- de_matrix[, cells_in, drop = FALSE]

  score <- rowMeans(x, na.rm = TRUE)

  if (standardize) {
    s <- matrixStats::rowSds(as.matrix(x), na.rm = TRUE)
    score <- score / (s + eps)
  }

  names(score) <- rownames(de_matrix)
  sort(score, decreasing = TRUE)
}

bic1_genes <- W[,5] > null_pos_W[5]
bic1_cells <- H[5,] > null_pos[5]
bic1_score <- score_bicluster(de_matrix, bic1_cells)


library(fgsea)
library(msigdbr)


scores <- score_bicluster(de_matrix, bic1_cells)

msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

res <- fgsea(pathways = pathways, stats = scores)
res <- res[order(res$padj), ]
head(res, 10)


scores_list <- lapply(1:ncol(W),
  function(i) {
    score_bicluster(de_matrix, (H[i,] > null_pos[i]))
})

res_list <- lapply(1:ncol(W),
  function(i) {
    fgsea(pathways = pathways, stats = scores_list[[i]])
})

all_pathways <- unique(unlist(lapply(res_list, function(x) x$pathway)))

nes_mat <- matrix(
  0,
  nrow = length(res_list),
  ncol = length(all_pathways),
  dimnames = list(names(res_list), all_pathways)
)

for (i in seq_along(res_list)) {
  res <- res_list[[i]]
  
  nes_mat[i, res$pathway] <- res$NES
}

library(pheatmap)

pheatmap(
  nes_mat,
  scale = "column",  # optional
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
d <- dist(nes_mat)
hc <- hclust(d)

plot(hc)

cor_mat <- cor(t(nes_mat))
pheatmap(cor_mat)




nes_mat_reduced <- nes_mat[, !grepl("INTERFERON", colnames(nes_mat))]
nes_var <- nes_mat_reduced[, apply(nes_mat_reduced, 2, var) > 0.2]
pheatmap::pheatmap(nes_var)

pca <- prcomp(nes_var, scale. = TRUE)
plot(pca$x[,1], pca$x[,2], pch = 19)
text(pca$x[,1], pca$x[,2], labels = rownames(nes_var), pos = 3)

plot(pca$x[,1], pca$x[,2], pch = 19,
     xlab = "PC1 (epithelial ↔ mesenchymal)",
     ylab = "PC2 (activity / metabolism)")

text(pca$x[,1], pca$x[,2],
     labels = paste0("B", 1:12), pos = 3)

