library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(Matrix)
library(lemur)
library(tidyverse)
library(Polychrome)
library(scico)
library(patchwork)

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

pval_threshold <- 1.

nei <- nei[order(nei[,"pval"]),]
nei_sig <- nei[nei[,"pval"] < pval_threshold,]
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

tibble(umap = umap) %>%
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

###############################
# FABIA!
jaccard <- function(A, B){
  num <- length(intersect(A, B))
  denom <- length(union(A, B))
  num/denom
}


library(fabia)

de_matrix <- as.matrix(assay(fit$training_data, "DE"))
#de_matrix <- tanh(0.5 * de_matrix)
de_matrix <- scale(de_matrix, center = TRUE, scale = TRUE)


data <- data.frame(value = as.vector(de_matrix))
ggplot(data, aes(x=value)) + 
  geom_histogram()

#de_matrix <- as.matrix(de_matrix)

res <- fabia(
  de_matrix,
  p = 5,        # number of factors (start large!)
  alpha = 0.5,  # sparsity (features)
  cyc = 1000      # iterations
)
# res is your fabia result
# res is your fabia result
ic <- res@avini   # information content per factor, length p
plot(sort(ic, decreasing = TRUE), type = "b",
     xlab = "Factor rank", ylab = "Information content",
     main = "Scree of fabia information content")


Z <- res@Z   # factors x cells
L <- res@L   # genes x factors


biclusters <- extractBic(res, thresZ = .5)

k <- 5
l <- 1
th_z <- 0.

scores <- res@Z[k, ]          # factors × samples → row k = sample scores for factor k

bic <- biclusters$bicopp[k,5]$biynn
bic2 <- biclusters$bic[2,5]$biypn

asc <- sort(L[which(abs(L[,k])>0.),k])
desc <- sort(L[which(abs(L[,k])>0.),k], decreasing=TRUE)

#sel_gene <- which(abs(L[,k])>0)[l]# is CXCL8


#L[sel_gene, k]

genes <- biclusters$bic[k,3]$bixn
genes <- names(sort(L[genes, k], decreasing=TRUE))
sel_gene <- genes[l]
#sel_gene <- "GATA2"


p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")[sel_gene,]) |>
  #arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()


p3 <- tibble(
  umap1     = umap[, 1],
  umap2     = umap[, 2],
  # Cleaned up the indicator logic
  indicator = ifelse(rownames(umap) %in% bic, "In List", "Other")
) |>
  ggplot(aes(x = umap1, y = umap2)) +
    # Use an anonymous function \(d) to tell subset() to use the piped tibble
    geom_point(
      data = \(d) subset(d, indicator == "Other"), 
      color = "grey85", 
      size = 0.4
    ) +
    geom_point(
      data = \(d) subset(d, indicator == "In List"), 
      color = "#0056B9",
      size = 0.6
    ) +
    coord_fixed() +
    theme_minimal()


ind <- which(nei$name == sel_gene)
neigh <- unlist(nei$neighborhood[ind])

p4 <- tibble(
  umap1     = umap[, 1],
  umap2     = umap[, 2],
  # Cleaned up the indicator logic
  indicator_old = ifelse(rownames(umap) %in% neigh, "In Neigh", "Other"),
  indicator_new = ifelse(rownames(umap) %in% bic, "In Bic", "Other"),
  indicator_int = ifelse(rownames(umap) %in% intersect(neigh, bic), "In Both", "Other")
) |>
  ggplot(aes(x = umap1, y = umap2)) +
    # Use an anonymous function \(d) to tell subset() to use the piped tibble
    geom_point(
      data = \(d) subset(d, indicator_old == "Other"), 
      color = "grey85", 
      size = 0.4
    ) +
    geom_point(
      data = \(d) subset(d, indicator_old == "In Neigh"), 
      color = "#327049",
      size = 0.6
    ) +
    geom_point(
      data = \(d) subset(d, indicator_new == "In Bic"), 
      color = "#0056B9",
      size = 0.6
    ) +
    geom_point(
      data = \(d) subset(d, indicator_int == "In Both"), 
      color = "#6c181f",
      size = 0.6
    ) +
    coord_fixed() +
    theme_minimal()

p1 + p3

p1 + p4
length(which(abs(L[,k])>0.))
L[sel_gene, k]

