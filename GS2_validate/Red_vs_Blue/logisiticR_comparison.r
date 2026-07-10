# library
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Define Gene Lists and Colour Mapping
# , "MALAT1"
red_genes <- c("TFG", "TNNI1") # shared
green_genes <- c("PRELID1", "SSBP2") # unique to original
blue_genes <- c("LDHA", "STARD10") # unique to refactored
my_genes <- c(red_genes, green_genes, blue_genes)

# Generate the interleaved target sequence
target_order <- as.vector(rbind(rev(red_genes), green_genes, blue_genes))

# Construct a reference data frame for gene categorisation
gene_colours <- data.frame(
  Gene = my_genes,
  ColourGroup = c(
    rep("Red", length(red_genes)),
    rep("Green", length(green_genes)),
    rep("Blue", length(blue_genes))
  )
)

# 2. Load OG Data
# sce_p1 and sg_allgenes
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/old_switching_genes_data.RData") 
output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/Red_vs_Blue/"

# 3. Process Expression Data and Add Noise
# Use assay() for standard slot access in SingleCellExperiment objects
sce1_GOI <- sce_p1[rownames(sce_p1) %in% my_genes, ]
sce1_GOI_df <- as.data.frame(t(assay(sce1_GOI, "expdata")))

# Define a seed to maintain reproducible noise generation
set.seed(42) 
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce1_GOI_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce1_GOI_df), 
  ncol = ncol(sce1_GOI_df)
)
sce1_GOI_df_noisy <- sce1_GOI_df + noise_matrix

# 5. Process Binarised Data for Logistic Regression
pseudotime_vec <- sce_p1$Pseudotime
gene_matrix <- t(as.matrix(assay(sce_p1, "binary")[my_genes, , drop = FALSE]))
gene_df <- as.data.frame(gene_matrix)
gene_df$Pseudotime <- pseudotime_vec

# Pivot and join with colour metadata
plot_data_long <- gene_df %>%
  pivot_longer(
    cols = all_of(my_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  left_join(gene_colours, by = "Gene") %>%
  mutate(FacetLabel = factor(Gene, levels = target_order))


# 2. Generate Faceted Logistic Fit Plot
png(paste0(output_dir, "OG-logistics.png"), width = 12, height = 6, units = "in", res = 300)

ggplot(plot_data_long, aes(x = Pseudotime, y = Expression)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.3, colour = "black") +
  stat_smooth(
    aes(colour = ColourGroup), 
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE
  ) +
  scale_colour_manual(
    values = c("Red" = "red", "Green" = "green", "Blue" = "blue")
  ) +
  facet_wrap(~ FacetLabel) +
  # Apply base theme first
  theme_minimal() + 
  labs(
    title = "Original implementation: Logistic Regression Fits",
    x = "Pseudotime",
    y = "Probability of Expression"
  ) +
  # Apply specific visual overrides second
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none" 
  )

dev.off()


# now show the logistic regression fits for the refactored implementation

# load the refactored data
# sce2 and gs2_allgenes
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/refactored_switching_genes_data.RData")

# 3. Process Expression Data and Add Noise
# Use assay() for standard slot access in SingleCellExperiment objects
sce2_GOI <- sce2[rownames(sce2) %in% my_genes, ]
sce2_GOI_df <- as.data.frame(t(assay(sce2_GOI, "expdata")))


# Define a seed to maintain reproducible noise generation
set.seed(42) 
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce2_GOI_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce2_GOI_df), 
  ncol = ncol(sce2_GOI_df)
)
sce2_GOI_df_noisy <- sce2_GOI_df + noise_matrix

# 5. Process Binarised Data for Logistic Regression
pseudotime_vec <- sce2$Pseudotime
gene_matrix <- t(as.matrix(assay(sce2_GOI, "binary")[my_genes, , drop = FALSE]))
gene_df <- as.data.frame(gene_matrix)
gene_df$Pseudotime <- pseudotime_vec

# Pivot and join with colour metadata
plot_data_long <- gene_df %>%
  pivot_longer(
    cols = all_of(my_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  left_join(gene_colours, by = "Gene") %>%
  mutate(FacetLabel = factor(Gene, levels = target_order))


# 2. Generate Faceted Logistic Fit Plot
png(paste0(output_dir, "GS2-logistics.png"), width = 12, height = 6, units = "in", res = 300)

ggplot(plot_data_long, aes(x = Pseudotime, y = Expression)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.3, colour = "black") +
  stat_smooth(
    aes(colour = ColourGroup), 
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE
  ) +
  scale_colour_manual(
    values = c("Red" = "red", "Green" = "green", "Blue" = "blue")
  ) +
  facet_wrap(~ FacetLabel) +
  # Apply base theme first
  theme_minimal() + 
  labs(
    title = "Refactored implementation: Logistic Regression Fits",
    x = "Pseudotime",
    y = "Probability of Expression"
  ) +
  # Apply specific visual overrides second
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none" 
  )

dev.off()