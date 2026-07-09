# library
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Define Gene Lists and Colour Mapping
# , "MALAT1"
red_genes <- c("TNNI1", "TFG") # shared
green_genes <- c("PRELID1", "SSBP2") # unique to original
blue_genes <- c("LDHA", "STARD10") # unique to refactored
my_genes <- c(red_genes, green_genes, blue_genes)

# Construct a reference data frame for gene categorisation
gene_colours <- data.frame(
  Gene = my_genes,
  ColourGroup = c(
    rep("Red", length(red_genes)),
    rep("Green", length(green_genes)),
    rep("Blue", length(blue_genes))
  )
)

output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/Red_vs_Blue/"

# 2. Load OG Data
# sce_p1 and sg_allgenes
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/old_switching_genes_data.RData") 

# 3. Process Expression Data and Add Noise
sce1_GOI <- sce_p1[rownames(sce_p1) %in% my_genes, ]
sce1_GOI_df <- as.data.frame(t(assay(sce1_GOI, "expdata")))

set.seed(42) 
# noise_matrix <- matrix(
#   rnorm(length(as.matrix(sce1_GOI_df)), mean = 0, sd = 0.1),
#   nrow = nrow(sce1_GOI_df), 
#   ncol = ncol(sce1_GOI_df)
# )
# sce1_GOI_df_noisy <- sce1_GOI_df + noise_matrix
sce1_GOI_df_noisy <- sce1_GOI_df 

# Pivot empirical data for faceting
expr_long <- sce1_GOI_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression") %>%
  left_join(gene_colours, by = "Gene")

# 4. Pre-calculate GMM Densities (Mixture calculation removed)
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

for (gene in my_genes) {
  params <- rowData(sce_p1)[gene, ]
  
  if (nrow(params) == 0) {
    warning(paste("Parameters for", gene, "not found. Skipping."))
    next
  }
  
  comp1_y <- params$lambda1 * dnorm(x_seq, mean = params$mu1, sd = params$sigma1)
  comp2_y <- params$lambda2 * dnorm(x_seq, mean = params$mu2, sd = params$sigma2)
  
  density_list[[gene]] <- data.frame(
    Gene = gene,
    x = x_seq,
    Comp1 = comp1_y,
    Comp2 = comp2_y,
    Root = params$root
  )
}

# 1. Construct the density data frame from the loop output
density_df <- bind_rows(density_list) %>%
  left_join(gene_colours, by = "Gene")

# Dynamically interleave the arrays: Red 1, Green 1, Blue 1, Red 2, Green 2, Blue 2
target_order <- as.vector(rbind(red_genes, green_genes, blue_genes))

# Enforce the order in both data frames
expr_long <- expr_long %>%
  mutate(Gene = factor(Gene, levels = target_order))

density_df <- density_df %>%
  mutate(Gene = factor(Gene, levels = target_order))

# 5. Generate Faceted Plot
png(paste0(output_dir, "OG_BINexpression_gmm.png"), width = 14, height = 10, units = "in", res = 300)

ggplot() +
  # Empirical Histogram Layer (Increased transparency, lighter borders)
  geom_histogram(
    data = expr_long,
    aes(x = Expression, y = after_stat(density), fill = ColourGroup),
    binwidth = 0.1, 
    colour = "darkgrey", 
    alpha = 0.3 
  ) +
  # Component 1 Layer
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp1), 
    linetype = "dashed", 
    colour = "red",
    linewidth = 0.8
  ) +
  # Component 2 Layer
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp2), 
    linetype = "dashed", 
    colour = "blue",
    linewidth = 0.8
  ) +
  # Decision Boundary Layer
  geom_vline(
    data = density_df, 
    aes(xintercept = Root), 
    linetype = "dotted", 
    linewidth = 1, 
    colour = "darkgreen"
  ) +
  # Structural Mapping
  facet_wrap(~ Gene, scales = "free_y") +
  # Updated to paler fill colours
  scale_fill_manual(values = c(
    "Red" = "lightcoral", 
    "Green" = "lightgreen", 
    "Blue" = "lightblue"
  )) +
  theme_minimal() +
  labs(
    title = "Original Implementation of Bimodal Gaussian Mixture Models",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

dev.off()


# repeat using the refactored data
# load the refactored data
# sce2 and sg2_allgenes
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/refactored_switching_genes_data.RData")
# 3. Process Expression Data and Add Noise
sce2_GOI <- sce2[rownames(sce2) %in% my_genes, ]
sce2_GOI_df <- as.data.frame(t(assay(sce2_GOI, "expdata")))

set.seed(42) 
# noise_matrix <- matrix(
#   rnorm(length(as.matrix(sce1_GOI_df)), mean = 0, sd = 0.1),
#   nrow = nrow(sce1_GOI_df), 
#   ncol = ncol(sce1_GOI_df)
# )
# sce1_GOI_df_noisy <- sce1_GOI_df + noise_matrix
sce2_GOI_df_noisy <- sce2_GOI_df 

# Pivot empirical data for faceting
expr_long <- sce2_GOI_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression") %>%
  left_join(gene_colours, by = "Gene")

# 4. Pre-calculate GMM Densities (Mixture calculation removed)
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

for (gene in my_genes) {
  params <- rowData(sce2)[gene, ]

  if (nrow(params) == 0) {
    warning(paste("Parameters for", gene, "not found. Skipping."))
    next
  }
  
  comp1_y <- params$lambda1 * dnorm(x_seq, mean = params$mu1, sd = params$sigma1)
  comp2_y <- params$lambda2 * dnorm(x_seq, mean = params$mu2, sd = params$sigma2)
  
  density_list[[gene]] <- data.frame(
    Gene = gene,
    x = x_seq,
    Comp1 = comp1_y,
    Comp2 = comp2_y,
    Root = params$root
  )
}

# 1. Construct the density data frame from the loop output
density_df <- bind_rows(density_list) %>%
  left_join(gene_colours, by = "Gene")

# Dynamically interleave the arrays: Red 1, Green 1, Blue 1, Red 2, Green 2, Blue 2
target_order <- as.vector(rbind(red_genes, green_genes, blue_genes))

# Enforce the order in both data frames
expr_long <- expr_long %>%
  mutate(Gene = factor(Gene, levels = target_order))

density_df <- density_df %>%
  mutate(Gene = factor(Gene, levels = target_order))

# 5. Generate Faceted Plot
png(paste0(output_dir, "GS2_BINexpression_gmm.png"), width = 14, height = 10, units = "in", res = 300)

ggplot() +
  # Empirical Histogram Layer (Increased transparency, lighter borders)
  geom_histogram(
    data = expr_long,
    aes(x = Expression, y = after_stat(density), fill = ColourGroup),
    binwidth = 0.1, 
    colour = "darkgrey", 
    alpha = 0.3 
  ) +
  # Component 1 Layer
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp1), 
    linetype = "dashed", 
    colour = "red",
    linewidth = 0.8
  ) +
  # Component 2 Layer
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp2), 
    linetype = "dashed", 
    colour = "blue",
    linewidth = 0.8
  ) +
  # Decision Boundary Layer
  geom_vline(
    data = density_df, 
    aes(xintercept = Root), 
    linetype = "dotted", 
    linewidth = 1, 
    colour = "darkgreen"
  ) +
  # Structural Mapping
  facet_wrap(~ Gene, scales = "free_y") +
  # Updated to paler fill colours
  scale_fill_manual(values = c(
    "Red" = "lightcoral", 
    "Green" = "lightgreen", 
    "Blue" = "lightblue"
  )) +
  theme_minimal() +
  labs(
    title = "Refactored Implementation of Bimodal Gaussian Mixture Models",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

dev.off()