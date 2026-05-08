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


folder_data <- "/home/herbermann/data/coculture/"
file_name <- "coculture_effect_in_AML_5000"
wd <- "/Users/herbermann/code/lemur_biclusters"
cache_dir <- file.path(wd, "cache")

# load or compute lemur fit
lemur_fit_file <- file.path(cache_dir, "coculture_lemur_fit.rds")
if (file.exists(lemur_fit_file)) {
  sce_all <- readRDS(paste0(folder_data, "sce_object_", file_name, ".rds"))
  lemur_fit <- readRDS(lemur_fit_file)
} else {
    
  sce_all <- readRDS(paste0(folder_data, "sce_object_", file_name, ".rds"))
    
  n_latents <- 20
  lemur_fit <- lemur(sce_all, ~coculture, n_embedding = n_latents, test_fraction = 0.5)
  lemur_fit <- align_by_grouping(lemur_fit, grouping = colData(fit)[, "clusters_renamed"])
  lemur_fit <- test_de(lemur_fit, contrast = cond(coculture = "Coculture") - cond(coculture = "BM_only"))
    
  saveRDS(lemur_fit, lemur_fit_file)
}


# Load or compute modified and split matrix
filename <- file.path(cache_dir, "coculture_cached_mat.rds")
if (file.exists(filename)) {
  mat <- readRDS(filename)
}else{
mat <- preprocess_assay_splits(
  lemur_fit,
  "DE",
  transform = "tanh",
  clip_max = "quantile",
  clip_max_value = 0.95
)}



filename <- file.path(cache_dir, "coculture_explore_factors.rds")
if (file.exists(filename)) {
  explore <- readRDS(filename)
}else{
  explore <- explore_factors_progress(
    mat$train,
    c(4,6,8,10,12,14,15,16,18,20),
    n_restarts = 10,
    tol = 0.05
  )
}

print(explore$plot)

# Factors we fit
n_factor = 15
factor_fit <- fit_nmf(mat$train, k=n_factor, n_restarts = 20, tol = 1e-3, seed = 42)









