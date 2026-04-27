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



## ============================================================
## Pseudobulk DE + Hallmark GSEA baseline for labeled cell types
## Paired design: individual + condition
## ============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(pheatmap)
  library(Matrix)
})

set.seed(42)

## -------------------------------
## User-adjustable column names
## -------------------------------
patient_col   <- "individual"
condition_col <- "coculture"         # e.g. BM_only / Coculture
celltype_col  <- "clusters_renamed"  # your annotated cell types

ref_condition <- "BM_only"
alt_condition <- "Coculture"

min_samples_per_condition <- 3       # per cell type, after pseudobulk
min_patients_per_celltype  <- 3       # optional extra filter

## -------------------------------
## Helper: check metadata columns
## -------------------------------
stopifnot(all(c(patient_col, condition_col, celltype_col) %in% colnames(colData(sce_all))))

## Force factors and reference level
sce_all[[patient_col]]   <- factor(sce_all[[patient_col]])
sce_all[[condition_col]] <- factor(sce_all[[condition_col]],
                                   levels = c(ref_condition, alt_condition))
sce_all[[celltype_col]]  <- factor(sce_all[[celltype_col]])

## -------------------------------
## Optional sanity checks
## -------------------------------
cat("Cells per condition x patient:\n")
print(table(sce_all[[condition_col]], sce_all[[patient_col]]))

cat("\nCells per cell type x condition:\n")
print(table(sce_all[[celltype_col]], sce_all[[condition_col]]))

## ------------------------------------------------------------
## 1) Pseudobulk aggregation by individual x cell type x condition
## ------------------------------------------------------------
pb <- aggregateAcrossCells(
  sce_all,
  ids = DataFrame(
    individual = sce_all[[patient_col]],
    condition   = sce_all[[condition_col]],
    cell_type   = sce_all[[celltype_col]]
  )
)

pb$individual <- factor(pb$individual)
pb$condition  <- factor(pb$condition, levels = c(ref_condition, alt_condition))
pb$cell_type  <- factor(pb$cell_type)

cat("\nPseudobulk samples:\n")
print(dim(pb))
print(table(pb$cell_type, pb$condition))

## ------------------------------------------------------------
## 2) Hallmark pathways
## ------------------------------------------------------------
msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

## ------------------------------------------------------------
## 3) Run paired pseudobulk DE for one cell type
## ------------------------------------------------------------
run_de_one_celltype <- function(pb, ct,
                                ref_condition = "BM_only",
                                alt_condition = "Coculture",
                                min_samples_per_condition = 3) {
  idx <- pb$cell_type == ct
  y <- pb[, idx]

  meta <- as.data.frame(colData(y))
  meta$individual <- factor(meta$individual)
  meta$condition  <- factor(meta$condition, levels = c(ref_condition, alt_condition))

  tab_cond <- table(meta$condition)
  if (any(tab_cond[c(ref_condition, alt_condition)] < min_samples_per_condition)) {
    return(NULL)
  }

  ## Need variation in condition and enough replicates
  if (nlevels(droplevels(meta$condition)) < 2) return(NULL)
  if (length(unique(meta$individual)) < min_patients_per_celltype) return(NULL)

  counts <- assay(y, "counts")

  dge <- DGEList(counts = counts)
  design <- model.matrix(~ individual + condition, data = meta)

  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  if (nrow(dge) < 50) return(NULL)  # too few genes left is usually not useful

  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)

  fit <- glmQLFit(dge, design)

  ## With ref_condition = BM_only and alt_condition = Coculture,
  ## the condition coefficient should be named "conditionCoculture"
  coef_name <- grep("^condition", colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    stop("Could not uniquely identify condition coefficient in design matrix.")
  }

  qlf <- glmQLFTest(fit, coef = coef_name)
  tt <- topTags(qlf, n = Inf)$table

  ## Signed ranking statistic for fgsea
  pval <- pmax(tt$PValue, .Machine$double.xmin)
  ranks <- sign(tt$logFC) * (-log10(pval))
  names(ranks) <- rownames(tt)
  ranks <- sort(ranks, decreasing = TRUE)

  ## fgsea
  gsea <- fgsea(
    pathways = pathways,
    stats    = ranks,
    minSize  = 10,
    maxSize  = 500,
    nperm    = 10000
  )

  gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]

  list(
    cell_type = ct,
    de_table  = tt,
    ranks     = ranks,
    gsea      = gsea,
    n_samples = ncol(y),
    n_patients = length(unique(meta$individual))
  )
}

## ------------------------------------------------------------
## 4) Run across cell types
## ------------------------------------------------------------
cell_types <- levels(pb$cell_type)

results <- lapply(cell_types, function(ct) {
  message("Running: ", ct)
  out <- tryCatch(
    run_de_one_celltype(
      pb,
      ct = ct,
      ref_condition = ref_condition,
      alt_condition = alt_condition,
      min_samples_per_condition = min_samples_per_condition
    ),
    error = function(e) {
      message("  failed: ", conditionMessage(e))
      NULL
    }
  )
  out
})
names(results) <- cell_types

results <- Filter(Negate(is.null), results)

cat("\nCell types successfully analyzed:\n")
print(names(results))

## ------------------------------------------------------------
## 5) Build NES matrix (cell type x pathway)
## ------------------------------------------------------------
all_pathways <- unique(unlist(lapply(results, function(x) x$gsea$pathway)))

nes_mat_cluster <- matrix(
  NA_real_,
  nrow = length(results),
  ncol = length(all_pathways),
  dimnames = list(names(results), all_pathways)
)

padj_mat_cluster <- matrix(
  NA_real_,
  nrow = length(results),
  ncol = length(all_pathways),
  dimnames = list(names(results), all_pathways)
)

for (i in seq_along(results)) {
  res <- results[[i]]
  g <- res$gsea
  nes_mat_cluster[i, g$pathway] <- g$NES
  padj_mat_cluster[i, g$pathway] <- g$padj
}

## ------------------------------------------------------------
## 6) Basic summary table
## ------------------------------------------------------------
summary_df <- do.call(
  rbind,
  lapply(results, function(x) {
    g <- x$gsea
    data.frame(
      cell_type = x$cell_type,
      n_samples = x$n_samples,
      n_patients = x$n_patients,
      n_sig_padj_0.05 = sum(g$padj < 0.05, na.rm = TRUE),
      top_pathway = if (nrow(g) > 0) g$pathway[1] else NA_character_,
      top_NES = if (nrow(g) > 0) g$NES[1] else NA_real_,
      top_padj = if (nrow(g) > 0) g$padj[1] else NA_real_,
      row.names = NULL
    )
  })
)

print(summary_df[order(summary_df$n_sig_padj_0.05, decreasing = TRUE), ])

sig_table <- do.call(
  rbind,
  lapply(results, function(x) {
    g <- x$gsea
    
    if (is.null(g) || nrow(g) == 0) return(NULL)
    
    g_sig <- g[g$padj < 0.05, ]
    if (nrow(g_sig) == 0) return(NULL)
    
    data.frame(
      cell_type = x$cell_type,
      n_samples = x$n_samples,
      n_patients = x$n_patients,
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

#saveRDS(sig_table, file = "rds/summary_classic_de.rds")

## ------------------------------------------------------------
## 7) Heatmap of NES values
## ------------------------------------------------------------
pheatmap(
  nes_mat_cluster,
  scale = "none",
  na_col = "grey90",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  border_color = NA,
  main = "Pseudobulk DE + Hallmark GSEA (cell-type baseline)"
)

## ------------------------------------------------------------
## 8) Optional: correlation among cell-type pathway profiles
## ------------------------------------------------------------
cor_mat_cluster <- cor(t(nes_mat_cluster), use = "pairwise.complete.obs")
pheatmap(
  cor_mat_cluster,
  scale = "none",
  na_col = "grey90",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  border_color = NA,
  main = "Correlation of Hallmark NES across cell types"
)

## ------------------------------------------------------------
## 9) Optional: compare to LEMUR bicluster NES matrix
##    (assumes you already have nes_mat_bicluster)
## ------------------------------------------------------------
# Common pathways only
# common_pw <- intersect(colnames(nes_mat_cluster), colnames(nes_mat_bicluster))
# cor_compare <- cor(
#   t(nes_mat_cluster[, common_pw, drop = FALSE]),
#   t(nes_mat_bicluster[, common_pw, drop = FALSE]),
#   use = "pairwise.complete.obs"
# )
# pheatmap(
#   cor_compare,
#   scale = "none",
#   border_color = NA,
#   main = "Cell-type baseline vs LEMUR bicluster NES correlation"
# )

## ------------------------------------------------------------
## 10) Save objects
## ------------------------------------------------------------
# saveRDS(results, "pseudobulk_gsea_results.rds")
# saveRDS(nes_mat_cluster, "nes_mat_cluster.rds")
# saveRDS(padj_mat_cluster, "padj_mat_cluster.rds")
# saveRDS(summary_df, "pseudobulk_gsea_summary.rds")
