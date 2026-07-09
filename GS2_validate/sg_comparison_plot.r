library(ggplot2)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)

#### comparison of the switching genes identified by both funcitons. ####
# seubset both to columns of interest
# [1] "geneID" "pseudoR2s"  "switch_at_timeidx"
og <- sg_allgenes[, c("geneID", "pseudoR2s", "switch_at_timeidx")]
gs2 <- sg2_allgenes[, c("geneID", "pseudoR2s", "switch_at_timeidx")]

#measure the overlap in geneID
matchingGenes <- intersect(og$geneID, gs2$geneID)
length(matchingGenes)
# 1317 the majority of genes are identified by both functions.

# subset to the genes which are unique to each function
uniqueGenes_OG <- setdiff(og$geneID, gs2$geneID)
# subset og to only include uniqueGenes
UQ_OG2 <- subset(og, geneID %in% uniqueGenes_OG)
# 24 genes unique to old
uniqueGenes_GS2 <- setdiff(gs2$geneID, og$geneID)
#subset gs2 to only inlcude uniqueGenes
UQ_GS2 <- subset(gs2, geneID %in% uniqueGenes_GS2)
# 25 genes unique to refactored.

# plot theswitching genes from both funcitons and highlight the unique ones from each funciton.

# merge sg_allgenes and sg2_allgenes
# update the feature_type column to say unique_original, unique_refactored or both.

UQ_sg2_allgenes <- subset(sg2_allgenes, geneID %in% uniqueGenes_GS2)

# Combine the full original dataset with unique genes from the refactored dataset
merged_allgenes <- rbind(sg_allgenes, UQ_sg2_allgenes)

# Classify the feature_type column using the predefined gene vectors
merged_allgenes$feature_type <- ifelse(
  merged_allgenes$geneID %in% uniqueGenes_OG, "unique_original",
  ifelse(merged_allgenes$geneID %in% uniqueGenes_GS2, "unique_refactored", "both")
)

# Convert feature_type to a factor for plotting
merged_allgenes$feature_type <- factor(
  merged_allgenes$feature_type,
  levels = c("both", "unique_original", "unique_refactored")
)

# 
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/merged_switching_genes_timeline.png", width = 10, height = 8, units = "in", res = 300)
plot_timeline_ggplot(merged_allgenes, timedata = sce2$Pseudotime, txtsize = 3)
dev.off()


#############


# 1. Define Gene Lists and Colour Mapping
red_genes <- c("TNNI1", "POU5F1", "MALAT1") # shared
green_genes <- c("PRELID1", "MARCKSL1") # unique to original
blue_genes <- c("LDHA", "PKM") # unique to refactored
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

# 2. Load Data
# sce_p1 and sg_allgenes
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/old_switching_genes_data.RData") 
output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/"

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
  mutate(FacetLabel = paste0(Gene, " (", ColourGroup, " List)"))

# 6. Generate Faceted Logistic Fit Plot
png(paste0(output_dir, "e-logistic_fit_my_genes.png"), width = 12, height = 6, units = "in", res = 300)

ggplot(plot_data_long, aes(x = Pseudotime, y = Expression)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.3, colour = "black") +
  stat_smooth(
    aes(colour = ColourGroup), 
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE
  ) +
  scale_colour_manual(
    values = c("Red" = "red", "Green" = "green", "Blue" = "blue")
  ) +
  facet_wrap(~ FacetLabel) +
  theme_minimal() +
  labs(
    title = "Binarised Expression vs Pseudotime",
    x = "Pseudotime",
    y = "Probability of Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none" 
  )

dev.off()


######### use expresion from SCE2 to see differnece in binarisation method. 

# Construct a reference data frame for gene categorisation
gene_colours <- data.frame(
  Gene = my_genes,
  ColourGroup = c(
    rep("Red", length(red_genes)),
    rep("Green", length(green_genes)),
    rep("Blue", length(blue_genes))
  )
)

# 2. Load the Refactored Data
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/refactored_switching_genes_data.RData")
output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/"

# 3. Process Expression Data from sce2 and Add Noise
sce2_GOI <- sce2[rownames(sce2) %in% my_genes, ]
sce2_GOI_df <- as.data.frame(t(assay(sce2_GOI, "expdata")))

set.seed(42) 
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce2_GOI_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce2_GOI_df), 
  ncol = ncol(sce2_GOI_df)
)
sce2_GOI_df_noisy <- sce2_GOI_df + noise_matrix

# 5. Process Binarised Data for Logistic Regression from sce2
pseudotime_vec <- sce2$Pseudotime
gene_matrix <- t(as.matrix(assay(sce2, "binary")[my_genes, , drop = FALSE]))
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
  mutate(FacetLabel = paste0(Gene, " (", ColourGroup, " List)"))

# 6. Generate Faceted Logistic Fit Plot (with gs2_ prefix)
# Add the gs2_ prefix to the file name
png(paste0(output_dir, "gs2_logistic_fit_my_genes.png"), width = 12, height = 6, units = "in", res = 300)

ggplot(plot_data_long, aes(x = Pseudotime, y = Expression)) +
  geom_jitter(height = 0.3, width = 0, alpha = 0.03, colour = "black") +
  stat_smooth(
    aes(colour = ColourGroup), 
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE
  ) +
  scale_colour_manual(
    values = c("Red" = "red", "Green" = "green", "Blue" = "blue")
  ) +
  facet_wrap(~ FacetLabel) +
  theme_minimal() +
  labs(
    title = "Binarised Expression vs Pseudotime (GS2)",
    x = "Pseudotime",
    y = "Probability of Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none" 
  )

dev.off()



##### plot pseudotime histogram

library(ggplot2)

# Construct the plotting data frame from the sce2 object
pseudotime_df <- data.frame(Pseudotime = sce2$Pseudotime)

# Open the PNG graphics device
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/gs2_pseudotime_histogram.png", width = 8, height = 6, units = "in", res = 300)

# Generate the histogram
ggplot(pseudotime_df, aes(x = Pseudotime)) +
  geom_histogram(
    binwidth = 2, 
    fill = "steelblue", 
    colour = "black"
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of Cells Across Pseudotime (GS2)",
    x = "Pseudotime",
    y = "Cell Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


######### repeat histogram with mixture models overlaid. - facceted.

gene_colours <- data.frame(
  Gene = my_genes,
  ColourGroup = c(
    rep("Red", length(red_genes)),
    rep("Green", length(green_genes)),
    rep("Blue", length(blue_genes))
  )
)

# 2. Load Data
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/old_switching_genes_data.RData")
output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/"

# 3. Process Expression Data and Add Noise
sce1_GOI <- sce_p1[rownames(sce_p1) %in% my_genes, ]
sce1_GOI_df <- as.data.frame(t(assay(sce1_GOI, "expdata")))

set.seed(42) 
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce1_GOI_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce1_GOI_df), 
  ncol = ncol(sce1_GOI_df)
)
sce1_GOI_df_noisy <- sce1_GOI_df + noise_matrix

# Pivot empirical data for faceting
expr_long <- sce1_GOI_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression") %>%
  left_join(gene_colours, by = "Gene")

# 4. Pre-calculate GMM Densities (Mixture calculation removed)
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

for (gene in my_genes) {
  params <- merged_allgenes[gene, ]
  
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

density_df <- bind_rows(density_list) %>%
  left_join(gene_colours, by = "Gene")

# 5. Generate Faceted Plot
png(paste0(output_dir, "faceted_expression_gmm.png"), width = 14, height = 10, units = "in", res = 300)

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
    title = "Faceted Gaussian Mixture Models",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

dev.off()


##### plot the distribution of roots. 
# as a histogram using ggplot2.
roots <- as.numeric(merged_allgenes$root)
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/merged_switching_genes_root_histogram.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data = data.frame(root = roots), aes(x = root)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Roots", x = "Root Value", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# which gene has the highest root value?
highest_root_gene <- merged_allgenes[which.max(merged_allgenes$root), ]
highest_root_gene


# see the full range of roots by using rowData(sce2)$root
roots2 <- rowData(sce2)$root
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/sce2_root_histogram.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data = data.frame(root = roots2), aes(x = root)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Roots in sce2", x = "Root Value", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# change the x lim to 0.1:6
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/sce2_root_histogram_zoomed.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data = data.frame(root = roots2), aes(x = root)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Roots in sce2 (Zoomed)", x = "Root Value", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0.3, 6) +
  ylim(0, 100)
dev.off()

# print the genes which have a root bigger than 0.4 and smaller than 1
medium_root_genes <- rowData(sce2)[rowData(sce2)$root >= 0.3 & rowData(sce2)$root < 1, ]
# randomly select 6 genes from medium_root_genes
set.seed(42)
random6_medium_root_genes <- sample(rownames(medium_root_genes), 6)
# plot the expression histograms with mixture models overlaid for these 6 genes

# 1. Filter and Sample Medium Root Genes
# Extract rowData as a data frame for easier manipulation
sce2_metadata <- as.data.frame(rowData(sce2))

# 2. Extract Data from sce2 and Add Noise
sce2_random6 <- sce2[rownames(sce2) %in% random6_medium_root_genes, ]
sce2_random6_df <- as.data.frame(t(assay(sce2_random6, "expdata")))

set.seed(42)
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce2_random6_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce2_random6_df),
  ncol = ncol(sce2_random6_df)
)
sce2_random6_df_noisy <- sce2_random6_df + noise_matrix

# Pivot empirical data for faceting
expr_long <- sce2_random6_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression")

# 3. Pre-calculate GMM Densities
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

for (gene in random6_medium_root_genes) {
  # Extract parameters directly from sce2 metadata
  params <- sce2_metadata[gene, ]
  
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

density_df <- bind_rows(density_list)

# 4. Generate the Faceted Plot
png(paste0(output_dir, "gs2_random6_medium_root_gmm.png"), width = 12, height = 8, units = "in", res = 300)

p <- ggplot() +
  # Empirical Histogram Layer
  geom_histogram(
    data = expr_long,
    aes(x = Expression, y = after_stat(density)),
    binwidth = 0.1, 
    fill = "lightgrey",
    colour = "darkgrey", 
    alpha = 0.4
  ) +
  # Component 1 Layer (Red Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp1), 
    linetype = "dotted", 
    colour = "red",
    linewidth = 1
  ) +
  # Component 2 Layer (Blue Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp2), 
    linetype = "dotted", 
    colour = "blue",
    linewidth = 1
  ) +
  # Decision Boundary Layer (Root)
  geom_vline(
    data = density_df, 
    aes(xintercept = Root), 
    linetype = "dashed", 
    linewidth = 1, 
    colour = "darkgreen"
  ) +
  # Create a 2x3 grid using facet_wrap
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(
    title = "Mixture Models for Randomly Sampled Medium-Root Genes (GS2)",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
dev.off()


### # print the genes which have a root bigger than 1 and smaller than Inf
# big_root_genes <- rowData(sce2)[rowData(sce2)$root >= 1 & rowData(sce2)$root <= 2.5, ]

# Filter for roots between 1 and 2.5
big_root_df <- sce2_metadata[
  !is.na(sce2_metadata$root) & 
  sce2_metadata$root >= 1 & 
  sce2_metadata$root <= 2.5, 
]

# Extract the target gene names as a character vector
big_root_genes <- rownames(big_root_df)

output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/"

# 2. Extract Data from sce2 and Add Noise
sce2_big_root_subset <- sce2[rownames(sce2) %in% big_root_genes, ]
sce2_big_root_df <- as.data.frame(t(assay(sce2_big_root_subset, "expdata")))

set.seed(42)
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce2_big_root_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce2_big_root_df),
  ncol = ncol(sce2_big_root_df)
)
sce2_big_root_df_noisy <- sce2_big_root_df + noise_matrix

# Pivot empirical data for faceting
expr_long <- sce2_big_root_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression")

# 3. Pre-calculate GMM Densities
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

# FIX: Iterate over the newly defined big_root_genes vector
for (gene in big_root_genes) {
  # Extract parameters directly from sce2 metadata
  params <- sce2_metadata[gene, ]
  
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

density_df <- bind_rows(density_list)

# 4. Generate the Faceted Plot
png(paste0(output_dir, "gs2_big_root_genes_gmm.png"), width = 14, height = 10, units = "in", res = 300)

p <- ggplot() +
  # Empirical Histogram Layer
  geom_histogram(
    data = expr_long,
    aes(x = Expression, y = after_stat(density)),
    binwidth = 0.1, 
    fill = "lightgrey",
    colour = "darkgrey", 
    alpha = 0.4
  ) +
  # Component 1 Layer (Red Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp1), 
    linetype = "dotted", 
    colour = "red",
    linewidth = 1
  ) +
  # Component 2 Layer (Blue Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp2), 
    linetype = "dotted", 
    colour = "blue",
    linewidth = 1
  ) +
  # Decision Boundary Layer (Root)
  geom_vline(
    data = density_df, 
    aes(xintercept = Root), 
    linetype = "dashed", 
    linewidth = 1, 
    colour = "darkgreen"
  ) +
  # Create a grid using facet_wrap
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  theme_minimal() +
  labs(
    title = "Mixture Models for Big-Root Genes (GS2)",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
dev.off()



####  se the range of switch_at_timeidx

switch_times <- merged_allgenes$switch_at_timeidx
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/merged_switching_genes_switch_time_histogram.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data = data.frame(switch_time = switch_times), aes(x = switch_time)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Switching Times", x = "Switching Time Index", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# looks good, and is whats shown in supplementary as better than switch de.




###### see the ditribution of means of the first distribution 
first_means <- sce2_metadata$mu1
# 
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/first_mean_histogram.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data = data.frame(first_mean = first_means), aes(x = first_mean)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of First Component Means", x = "Mean of First Component (mu1)", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off() 

png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/first_mean_histogram_zoomed.png", width = 8, height = 6, units = "in", res = 300)
ggplot(data = data.frame(first_mean = first_means), aes(x = first_mean)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of First Component Means (Zoomed)", x = "Mean of First Component (mu1)", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0.1, 5) +
  ylim(0, 500)
dev.off()

# 

### # print the genes which have a first mean bigger than 0.1 

# Filter for mu above 0.1
big_mu1_df <- sce2_metadata[
  !is.na(sce2_metadata$mu1) & 
  sce2_metadata$mu1 > 0.1, 
]

# Extract the target gene names as a character vector
big_mu1_genes <- rownames(big_mu1_df)

# Randomly sample exactly 9 genes
set.seed(42)
sampled_mu1_genes <- sample(big_mu1_genes, 9)

# 2. Extract Data from sce2 and Add Noise
# Isolate the data using the sampled 9 genes
sce2_sampled_mu1_subset <- sce2[rownames(sce2) %in% sampled_mu1_genes, ]
sce2_sampled_mu1_df <- as.data.frame(t(assay(sce2_sampled_mu1_subset, "expdata")))

set.seed(42)
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce2_sampled_mu1_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce2_sampled_mu1_df),
  ncol = ncol(sce2_sampled_mu1_df)
)
sce2_sampled_mu1_df_noisy <- sce2_sampled_mu1_df + noise_matrix

# Pivot empirical data for faceting
expr_long <- sce2_sampled_mu1_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression")

# 3. Pre-calculate GMM Densities
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

# Iterate strictly over the 9 sampled genes
for (gene in sampled_mu1_genes) {
  params <- sce2_metadata[gene, ]
  
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

density_df <- bind_rows(density_list)

# 4. Generate the Faceted Plot
png(paste0(output_dir, "gs2_sampled_mu1_genes_gmm.png"), width = 12, height = 10, units = "in", res = 300)

p <- ggplot() +
  # Empirical Histogram Layer
  geom_histogram(
    data = expr_long,
    aes(x = Expression, y = after_stat(density)),
    binwidth = 0.1, 
    fill = "lightgrey",
    colour = "darkgrey", 
    alpha = 0.4
  ) +
  # Component 1 Layer (Red Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp1), 
    linetype = "dotted", 
    colour = "red",
    linewidth = 1
  ) +
  # Component 2 Layer (Blue Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp2), 
    linetype = "dotted", 
    colour = "blue",
    linewidth = 1
  ) +
  # Decision Boundary Layer (Root)
  geom_vline(
    data = density_df, 
    aes(xintercept = Root), 
    linetype = "dashed", 
    linewidth = 1, 
    colour = "darkgreen"
  ) +
  # Create a 3x3 grid using facet_wrap
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(
    title = "Mixture Models for Sampled High-Mu1 Genes (GS2)",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
dev.off()


######

# Filter for mu1 strictly between 0.1 and 0.5
medium_mu1_df <- sce2_metadata[
  !is.na(sce2_metadata$mu1) & 
  sce2_metadata$mu1 > 0.1 & 
  sce2_metadata$mu1 < 0.7, 
]

# Extract the target gene names
medium_mu1_genes <- rownames(medium_mu1_df)

# Protect against datasets with fewer than 9 valid genes
n_sample <- min(9, length(medium_mu1_genes))

set.seed(42)
sampled_medium_mu1_genes <- sample(medium_mu1_genes, n_sample)

output_dir <- "/home/mtn1n22/scratch/thesisFigures/GS2_validate/"

# 2. Extract Data from sce2 and Add Noise
sce2_medium_mu1_subset <- sce2[rownames(sce2) %in% sampled_medium_mu1_genes, ]
sce2_medium_mu1_df <- as.data.frame(t(assay(sce2_medium_mu1_subset, "expdata")))

set.seed(42)
noise_matrix <- matrix(
  rnorm(length(as.matrix(sce2_medium_mu1_df)), mean = 0, sd = 0.1),
  nrow = nrow(sce2_medium_mu1_df),
  ncol = ncol(sce2_medium_mu1_df)
)
sce2_medium_mu1_df_noisy <- sce2_medium_mu1_df + noise_matrix

# Pivot empirical data for faceting
expr_long <- sce2_medium_mu1_df_noisy %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression")

# 3. Pre-calculate GMM Densities
x_seq <- seq(min(expr_long$Expression), max(expr_long$Expression), length.out = 500)
density_list <- list()

for (gene in sampled_medium_mu1_genes) {
  params <- sce2_metadata[gene, ]
  
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

density_df <- bind_rows(density_list)

# 4. Generate the Faceted Plot
png(paste0(output_dir, "gs2_sampled_medium_mu1_genes_gmm.png"), width = 12, height = 10, units = "in", res = 300)

p <- ggplot() +
  # Empirical Histogram Layer
  geom_histogram(
    data = expr_long,
    aes(x = Expression, y = after_stat(density)),
    binwidth = 0.1, 
    fill = "lightgrey",
    colour = "darkgrey", 
    alpha = 0.4
  ) +
  # Component 1 Layer (Red Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp1), 
    linetype = "dotted", 
    colour = "red",
    linewidth = 1
  ) +
  # Component 2 Layer (Blue Dotted)
  geom_line(
    data = density_df, 
    aes(x = x, y = Comp2), 
    linetype = "dotted", 
    colour = "blue",
    linewidth = 1
  ) +
  # Decision Boundary Layer (Root)
  geom_vline(
    data = density_df, 
    aes(xintercept = Root), 
    linetype = "dashed", 
    linewidth = 1, 
    colour = "darkgreen"
  ) +
  # Create a grid using facet_wrap
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(
    title = "Mixture Models for Medium-Mu1 Genes (GS2)",
    x = "Log-normalised Expression",
    y = "Density"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
dev.off()