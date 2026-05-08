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
set.seed(42)


folder_data <- "/Users/herbermann/code/data/coculture/"
file_name <- "coculture_effect_in_AML_5000"
wd <- "/Users/herbermann/code/lemur_biclusters"
cache_dir <- file.path(wd, "cache")

lemur_fit_file <- file.path(cache_dir, "lemur_essi_fit.rds")

FLAG_fit_lemur <- TRUE

if (FLAG_fit_lemur == TRUE){
  if (file.exists(lemur_fit_file)) {
    sce_all <- readRDS(paste0(folder_data, "sce_object_", file_name, ".RDS"))
    lemur_fit <- readRDS(lemur_fit_file)
  } else {
    
    sce_all <- readRDS(paste0(folder_data, "sce_object_", file_name, ".RDS"))
    
    n_latents <- 20
    lemur_fit <- lemur(sce_all, ~coculture, n_embedding = n_latents, test_fraction = 0.5)
    lemur_fit <- align_by_grouping(lemur_fit, grouping = colData(fit)[, "clusters_renamed"])
    lemur_fit <- test_de(lemur_fit, contrast = cond(coculture = "Coculture") - cond(coculture = "BM_only"))
    
    saveRDS(lemur_fit, lemur_fit_file)
  }
}


mat <- preprocess_assay_splits(
  lemur_fit,
  "DE",
  transform = "tanh",
  clip_max = "quantile",
  clip_max_value = 0.95
)