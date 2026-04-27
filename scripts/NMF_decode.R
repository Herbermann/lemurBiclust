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


folder_data <- "/Users/herbermann/code/data/lemur4decode"
filename <- "sce_Notch_clean.rds"

sce <- readRDS(file.path(folder_data, filename))


set.seed(42)

n_latents <- 20

fit = lemur(sce, design = ~ Batch + Perturbation, n_embedding = n_latents, test_fraction = 0.2) 
fit = align_harmony(fit) 

#| label: umap
umap = uwot::umap(t(fit$embedding), n_neighbors = 30, scale = "none")
reducedDim(fit, "UMAP_lemur") <- umap

#| label: umapdf
df = as_tibble(fit$colData) |> mutate(umap = umap)


#| label: umapcolors
celltype_colors =  c(
  "ISC" = "#FFE0F2",
  "EB" = "#FF78D5",
  "daEC" = "#81E6D4",
  "aEC" = "#1c9f78",
  "dpEC" = "#78CEFF",
  "pEC" = "#218BE2",
  "EEP" = "#fcdb03",
  "AstC-EE" = "#feb164",
  "Tk-EE" = "#fd9997",
  "EE-uncert" = "#e51a1a",
  "mEC" = "#BA90DE",
  "CC" = "#6a3c9b",
  "LFC" = "#C9F746",
  "NAAT1-MtnA enr" = "#853a38",
  "MtnA-Smvt enr" = "#FB3C7A",
  "unk" = "#808080"
)
#assert_that(is.factor(df$cell_type), setequal(levels(df$cell_type), names(celltype_colors)))

ggplot(df, aes(x = umap[,1], y = umap[,2])) +
  geom_point(aes(color = cell_type), size = 1, alpha = 0.5) +
  scale_color_manual(values = celltype_colors, na.value = celltype_colors["unk"]) + # unexpectedly, there are unannotated cells. Need to look into that.
  coord_fixed() +
  labs(color = "cell type") +
  theme(
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key.size = unit(0.6, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(vars(Batch))


fit_de = test_de(fit, contrast = cond(Perturbation = "N") - cond(Perturbation = "control"))

#### BICLUST GO!

############################################
library(Matrix)
# Development version from GitHub (requires Rcpp and RcppEigen)
library("RcppML")
packageVersion("RcppML")


select_hvg <- function(mat, n_top = 3000) {
  stopifnot(inherits(mat, "dgCMatrix"))
  
  # mean per gene
  gene_means <- Matrix::rowMeans(mat)
  
  # variance per gene (sparse-safe)
  gene_vars <- Matrix::rowMeans(mat^2) - gene_means^2
  
  # avoid division by zero
  gene_means[gene_means == 0] <- 1e-8
  
  # dispersion (variance normalized by mean)
  disp <- gene_vars / gene_means
  
  # rank genes
  top_idx <- order(disp, decreasing = TRUE)[seq_len(min(n_top, length(disp)))]
  
  mat[top_idx, ]
}


de_matrix <- Matrix(assay(fit_de$training_data, "DE"), sparse=TRUE)


de_matrix <- select_hvg(de_matrix, n_top = 3000)


de_pos <- de_matrix
de_pos@x <- pmax(de_pos@x, 0)

de_neg <- de_matrix
de_neg@x <- pmax(-de_neg@x, 0)

de_pos <- drop0(de_pos)
de_neg <- drop0(de_neg)

quant_pos <- quantile(de_pos@x, c(0.5, 0.75, 0.9, 0.95))
quant_neg <- quantile(de_neg@x, c(0.5, 0.75, 0.9, 0.95))

de_pos@x <- tanh(2 / quant_pos[4] * de_pos@x)
de_neg@x <- tanh(2 / quant_neg[4] * de_neg@x)

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

      fit <- RcppML::nmf(train_matrix, k = k, seed = i, L1 = L1, tol = tol)

      err <- if (use_cv) {
        RcppML::evaluate(fit, data = matrix, mask = test_mask, missing_only = TRUE)
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

    #print(p_error + p_stab)
  }

  list(results = results, plot = (p_error + p_stab))
}

with_progress({
  test <- find_factors(de, c(3, 4, 5, 6, 12, 13, 14, 15, 16), n_restarts=3, test_fraction = 0.1 , show_plot = TRUE, L1 = c(0.0 , 0.0) , tol = 1e-2)
})

test$plot

# we can try k = 15 here.

fit_nmf <- function(matrix,
              k,
              L1 = c(0, 0),
              n_restarts = 1,
              n_perm = 10,
              thr_h = 0.95,
              thr_w_pos = 0.95,
              thr_w_neg = 0.95
              )
  {

  biclusters <- nmf(matrix, k = k, seed = 1:n_restarts)
  
  n_genes <- nrow(matrix) / 2

  nulls <- lapply(seq_len(n_perm), function(i) {
    mat_perm <- t(apply(matrix, 1, sample))
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-2)
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

  list(biclusters = biclusters, null_h = null_threshold_h, null_w_pos = null_threshold_w_pos, null_w_neg = null_threshold_w_neg)
}


nmf_fit <- fit_nmf(de, k = 14, n_restart = 1, n_perm = 5, thr_h = c(0.9, 0.95, 0.99, 0.995, 0.999), thr_w_pos = c(0.9, 0.95, 0.99, 0.995, 0.999), thr_w_neg = c(0.9, 0.95, 0.99, 0.995, 0.999))

biclusters <- nmf_fit$biclusters
null_h <- nmf_fit$null_h
null_w_pos <- nmf_fit$null_w_pos
null_w_neg <- nmf_fit$null_w_neg

W <- biclusters@w
H <- biclusters@h

n_genes <- nrow(de)/2

W_pos <- W[1:n_genes, ]
W_neg  <- W[(n_genes+1):nrow(de), ]

rownames(W_pos) <- pos_orig
rownames(W_neg) <- neg_orig

null_h <- null_h[3,]
null_w_pos <- null_w_pos[3,]
null_w_neg <- null_w_neg[3,]


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

bic_genes <- Map(function(up, down) up | down, bic_up, bic_down)

summary(lapply(bic_genes, function(x) {which(x) == TRUE}))

k <- 2
l <- 1

de_expression <- de_matrix

genes_up <- sort(W_pos[,k], decreasing = TRUE)

genes_down <- sort(W_neg[,k], decreasing = TRUE)

sel_gene <- names(genes_up[l])
#sel_gene <- "GATA2"

umap <- uwot::umap(t(fit_de$training_data$embedding))
#umap <-uwot::umap(t(assay(fit$training_data,"DE")))

p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit_de$training_data, "DE")[sel_gene,]) |>
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

p1

p2



n_cells <- length(bic_cells[[1]])
n_bics  <- length(bic_cells)

cell_label <- sapply(seq_len(n_cells), function(i) {
  bics_here <- which(vapply(seq_len(n_bics), function(b) bic_cells[[b]][i], logical(1)))
  if (length(bics_here) == 0) NA_character_
  else paste0("BIC", bics_here[1])  # first-wins for overlapping cells
})

p3 <- tibble(x = umap[,1], y = umap[,2]) |>
  mutate(bic = factor(cell_label)) |>
  arrange(!is.na(bic)) |>          # NA (grey) points drawn first
  ggplot(aes(x, y)) +
    geom_point(aes(color = bic), size = 0.5) +
    scale_color_discrete(na.value = "grey85") +
    coord_fixed()

p3


#############
library(pheatmap)

jaccard_with_logic <- function(A, B) {
  stopifnot(length(A) == length(B))
  denom <- sum(A | B)
  if (denom == 0) return(NA_real_)
  sum(A & B) / denom
}

n_bic <- length(bic_cells)

jaccard_cells <- matrix(
  NA_real_,
  nrow = n_bic,
  ncol = n_bic,
  dimnames = list(paste0("Bic_", seq_len(n_bic)),
                  paste0("Bic_", seq_len(n_bic)))
)

for (i in seq_len(n_bic)) {
  for (j in seq_len(n_bic)) {
    jaccard_cells[i, j] <- jaccard_with_logic(bic_cells[[i]], bic_cells[[j]])
  }
}

n_bic <- length(bic_genes)

jaccard_genes <- matrix(
  NA_real_,
  nrow = n_bic,
  ncol = n_bic,
  dimnames = list(paste0("BIC", seq_len(n_bic)),
                  paste0("BIC", seq_len(n_bic)))
)

for (i in seq_len(n_bic)) {
  for (j in seq_len(n_bic)) {
    jaccard_genes[i, j] <- jaccard_with_logic(bic_genes[[i]], bic_genes[[j]])
  }
}

pheatmap(
  jaccard_cells,
  scale = "none",
#  cluster_rows = FALSE,
#  cluster_cols = FALSE
)
pheatmap(
  jaccard_genes,
  scale = "none",
#  cluster_rows = FALSE,
#  cluster_cols = FALSE
)

######

celltype_map <- setNames(sce$cell_type, colnames(sce))

#Investigate connection of biclusters to known annotations
# --- align all objects to the same cells ---
common_cells <- Reduce(intersect, c(
  list(names(celltype_map)),
  lapply(bic_cells, names)
))

celltype_map2 <- celltype_map[common_cells]
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



heat_df <- plot_df %>%
  select(celltype, bic, prop) %>%
  pivot_wider(names_from = bic, values_from = prop, values_fill = 0)
heat_mat <- as.matrix(heat_df[,-1])
rownames(heat_mat) <- heat_df$celltype

cols <- colorRampPalette(c("white", "#C6DBEF", "#6BAED6", "#2171B5"))(100)

pheatmap(
  heat_mat,
  color = cols,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Cell type composition across biclusters"
)



#########
library(fgsea)
library(msigdbr)

msig <- msigdbr(species = "Drosophila melanogaster", category = "C5", subcategory = "GO:BP")
pathways <- split(msig$gene_symbol, msig$gs_name)

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



scores_list <- lapply(1:ncol(W_pos),
  function(i) {
    score_bicluster(de_matrix, (H[i,] > null_h[i]))
})

res_list <- lapply(1:ncol(W_pos),
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

res_list_filtered <- lapply(seq_along(scores_list), function(i) {
  res <- fgsea(pathways = pathways, stats = scores_list[[i]])
  
  res <- res[res$padj < 0.05 & abs(res$NES) > 1.]
  res <- res[order(-abs(res$NES)), ]
  
  if (nrow(res) == 0) return(NULL)
  
  collapsed <- collapsePathways(res, pathways, stats = scores_list[[i]])
  res <- res[res$pathway %in% collapsed$mainPathways, ]
  
  res
})

all_pathways <- unique(unlist(lapply(res_list_filtered, function(x) {
  if (is.null(x)) return(NULL)
  x$pathway
})))
bic_names <- paste0("BIC", seq_along(res_list_filtered))

nes_mat <- matrix(
  0.,#NA_real_,
  nrow = length(res_list_filtered),
  ncol = length(all_pathways),
  dimnames = list(bic_names, all_pathways)
)

for (i in seq_along(res_list_filtered)) {
  res <- res_list_filtered[[i]]
  if (is.null(res)) next
  
  nes_mat[i, res$pathway] <- res$NES
}


shorten_names <- function(x, max_len = 40) {
  ifelse(nchar(x) > max_len,
         paste0(substr(x, 1, max_len), "..."),
         x)
}
colnames(nes_mat) <- shorten_names(colnames(nes_mat), 40)


pheatmap(
  nes_mat,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)


cor_mat <- cor(t(nes_mat), use = "pairwise.complete.obs")

pheatmap(
  cor_mat,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

#####


pca_nes <- prcomp(nes_mat, center = TRUE, scale. = TRUE)

pca_nes$x |>
  as.data.frame() |>
  rownames_to_column("bicluster") |>
  ggplot(aes(PC1, PC2, color = bicluster, label = bicluster)) +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.3, show.legend = FALSE) +
  theme_minimal()

pca_nes$x |>
  as.data.frame() |>
  rownames_to_column("bicluster") |>
  ggplot(aes(PC2, PC3, color = bicluster, label = bicluster)) +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.3, show.legend = FALSE) +
  theme_minimal()

summary(pca_nes)


pca_nes$rotation |>
  pheatmap(
    scale = 'none',
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    angle_col = 45,
    fontsize_row = 8
  )
