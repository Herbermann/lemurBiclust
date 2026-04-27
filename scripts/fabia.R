devtools::load_all("/Users/herbermann/code/lemur_biclusters/lemur")
library("lemur")
library("tidyverse")
library("SingleCellExperiment")
library("patchwork")

set.seed(42)

data("glioblastoma_example_data", package = "lemur")
glioblastoma_example_data

fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition, 
             n_embedding = 15, test_fraction = 0.5)

fit <- align_harmony(fit)

umap <- uwot::umap(t(fit$embedding))

fit <- compute_contrasts(fit, contrast = cond(condition = "panobinostat") - cond(condition = "ctrl"))

class(assay(fit, "DE"))

sel_gene <- "ENSG00000169429" # is CXCL8

p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")[sel_gene,]) |>
  arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()


neighborhoods <- find_de_neighborhoods(fit, group_by = vars(patient_id, condition))

as_tibble(neighborhoods) |>
left_join(as_tibble(rowData(fit)[,1:2]), by = c("name" = "gene_id")) |>
  relocate(symbol, .before = "name") |>
  arrange(pval) |>
  head(5)

####################
library(fabia)

jaccard <- function(A, B){
  num <- length(intersect(A, B))
  denom <- length(union(A, B))
  num/denom
}


w <- 0.5

b <- 10
a <- 0.4

smoother <- function(x, b, a){
  1/(1 + exp(-b * (abs(x)-a)))
}


de_matrix <- as.matrix(assay(fit, "DE"))
scale(de_matrix, center=TRUE, scale=TRUE)
#de_matrix <- tanh(w * de_matrix)

# Allow for sigmoidal noise suppression
#de_matrix <- de_matrix * smoother(de_matrix, b, a)

#de_matrix <- scale(de_matrix)

library(matrixStats)
#de_matrix <- de_matrix - rowMedians(de_matrix)

# 2. Scale
#de_matrix <- de_matrix / rowSds(de_matrix)


data <- data.frame(value = as.vector(de_matrix))
ggplot(data, aes(x=value)) + 
  geom_histogram()

#de_matrix <- as.matrix(de_matrix)

res <- fabia(
  de_matrix,
  p = 3,        # number of factors (start large!)
  alpha = 0.5,  # sparsity (features)
  cyc = 500      # iterations
)

# res is your fabia result
ic <- res@avini   # information content per factor, length p
plot(sort(ic, decreasing = TRUE), type = "b",
     xlab = "Factor rank", ylab = "Information content",
     main = "Scree of fabia information content")




Z <- res@Z   # factors x cells
L <- res@L   # genes x factors


biclusters <- extractBic(res, thresZ = .5)

k <- 1
l <- 1

bic <- biclusters$bic[k,5]$biypn
bicopp <- biclusters$bicopp[k,5]$biynn
#asc <- sort(L[which(abs(L[,k])>0.),k])
#desc <- sort(L[which(abs(L[,k])>0.),k], decreasing=TRUE)

#sel_gene <- which(abs(L[,k])>0)[l]# is CXCL8


#L[sel_gene, k]

genes <- biclusters$bic[k,3]$bixn
genes_opp <- biclusters$bicopp[k,3]$bixn

genes2 <- biclusters$bic[2, 3]$bixn

#genes <- names(sort(L[genes, k], decreasing=TRUE))
loadings_subset <- res@L[genes, ]

real_k <- which.max(colMeans(abs(loadings_subset)))

sel_gene <- genes[l]

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

ind <- which(neighborhoods$name == sel_gene)
neigh <- unlist(neighborhoods$neighborhood[ind])

p4 <- tibble(
  umap1     = umap[, 1],
  umap2     = umap[, 2],
  # Cleaned up the indicator logic
  indicator_old = ifelse(rownames(umap) %in% neigh, "In Neigh", "Other"),
  indicator_new = ifelse(rownames(umap) %in% bic, "In Bic", "Other"),
  indicator_opp = ifelse(rownames(umap) %in% bicopp, "In Bicopp", "Other"),
  indicator_int = ifelse(rownames(umap) %in% intersect(neigh, bic), "In Both", "Other"),
  indicator_intopp = ifelse(rownames(umap) %in% intersect(neigh, bicopp), "In Both opp", "Other")
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
      data = \(d) subset(d, indicator_opp == "In Bicopp"), 
      color = "#2ab2c4",
      size = 0.6
    ) +
  
    geom_point(
      data = \(d) subset(d, indicator_int == "In Both"), 
      color = "#6c181f",
      size = 0.6
    ) +
  
    geom_point(
      data = \(d) subset(d, indicator_intopp == "In Both opp"), 
      color = "#e89300",
      size = 0.6
    ) +
  
    coord_fixed() +
    theme_minimal()


correct_k <- k


cell_loadings <- res@Z[correct_k, ]
cell_loadings_names <- names(cell_loadings)



# 1. Get the raw Z vector for your factor
# We use drop=TRUE to ensure it's a simple numeric vector
z_raw <- res@Z[correct_k, ]

# 2. Assign names from the ORIGINAL matrix used in fabia
# FABIA doesn't always keep names in the @Z slot, so we pull from the input
names(z_raw) <- colnames(de_matrix)

# 3. Align to the UMAP rownames
# This handles the 5000 vs 2500 mismatch by filling missing cells with NA
z_aligned <- z_raw[rownames(umap)]

# 4. Build the plot
p_check <- tibble(
  u1 = umap[, 1],
  u2 = umap[, 2],
  z  = z_aligned,
  cell_id = rownames(umap)
)

ggplot(p_check, aes(x = u1, y = u2)) +
  # Background: All cells in UMAP
  geom_point(color = "grey85", size = 0.4) +
  # Foreground: Only cells that were actually in the FABIA run
  geom_point(
    data = \(d) subset(d, !is.na(z)),
    aes(color = z), 
    size = 0.6
  ) +
  scale_color_gradient2(low = "#FFD800", mid = "white", high = "#0056B9", midpoint = 0) +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Factor", correct_k, "Continuous Loading"))



p1 + p3

p1 + p4

L[sel_gene, correct_k]




figure_wrapper <- function(biclusters, k, neighborhoods){

  bic <- biclusters$bic[k,5]$biypn
  bicopp <- biclusters$bicopp[k,5]$biynn
  genes <- biclusters$bic[k,3]$bixn
  genes_opp <- biclusters$bicopp[k,3]$bixn


  for (l in 1:length(genes)) {
    
    sel_gene <- genes[l]

    p<- tibble(umap = umap) |>
      mutate(de = assay(fit, "DE")[sel_gene,]) |>
      #arrange(abs(de)) |>
      ggplot(aes(x = umap[,1], y = umap[,2])) +
        geom_point(aes(color = de), size = 0.5) +
        scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
        coord_fixed()


    ind <- which(neighborhoods$name == sel_gene)
    neigh <- unlist(neighborhoods$neighborhood[ind])

    pp <- tibble(
      umap1     = umap[, 1],
      umap2     = umap[, 2],
      # Cleaned up the indicator logic
      indicator_old = ifelse(rownames(umap) %in% neigh, "In Neigh", "Other"),
      indicator_new = ifelse(rownames(umap) %in% bic, "In Bic", "Other"),
      indicator_opp = ifelse(rownames(umap) %in% bicopp, "In Bicopp", "Other"),
      indicator_int = ifelse(rownames(umap) %in% intersect(neigh, bic), "In Both", "Other"),
      indicator_intopp = ifelse(rownames(umap) %in% intersect(neigh, bicopp), "In Both opp", "Other")
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
          color = "#e89300",
          size = 0.6
        ) +    
        coord_fixed() +
        theme_minimal()

    ppp <- tibble(
      umap1     = umap[, 1],
      umap2     = umap[, 2],
      # Cleaned up the indicator logic
      indicator_old = ifelse(rownames(umap) %in% neigh, "In Neigh", "Other"),
      indicator_new = ifelse(rownames(umap) %in% bic, "In Bic", "Other"),
      indicator_opp = ifelse(rownames(umap) %in% bicopp, "In Bicopp", "Other"),
      indicator_int = ifelse(rownames(umap) %in% intersect(neigh, bic), "In Both", "Other"),
      indicator_intopp = ifelse(rownames(umap) %in% intersect(neigh, bicopp), "In Both opp", "Other")
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
          data = \(d) subset(d, indicator_opp == "In Bicopp"), 
          color = "#0056B9",
          size = 0.6
        ) +
        geom_point(
          data = \(d) subset(d, indicator_intopp == "In Both opp"), 
          color = "#e89300",
          size = 0.6
        ) +
      
        coord_fixed() +
        theme_minimal()
    
        print(p + pp + ppp)

  }

}


jac_them <- function(biclusters, k, neighborhoods) {

  bic <- biclusters$bic[k,5]$biypn
  bicopp <- biclusters$bicopp[k,5]$biynn
  genes <- biclusters$bic[k,3]$bixn
  genes_opp <- biclusters$bicopp[k,3]$bixn

  df <- tibble(gene = genes) %>%
    mutate(
      # Find index for each gene
      match_ind = map_int(gene, ~ {
        idx <- which(neighborhoods$name == .x)
        if(length(idx) > 0) idx[1] else NA_integer_
      }),
      # Extract neighborhood and calculate Jaccard
      jac_bic = map_dbl(match_ind, ~ {
        if(is.na(.x)) return(NA)
        jaccard(bic, unlist(neighborhoods$neighborhood[.x]))
      }),
      jac_opp = map_dbl(match_ind, ~ {
        if(is.na(.x)) return(NA)
        jaccard(bicopp, unlist(neighborhoods$neighborhood[.x]))
      })
    )
  df
}



