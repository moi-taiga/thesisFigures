### compare the Binarisation funcitons outputs ####
# using kappa 

# assays(sce)$expdata is not sparse (dgCMatrix).
# Coercing expression data to dgCMatrix for efficiency...
# 305 genes are expressed in 5 or fewer cells and will be binarized to 0's.
# Fitting mixture models for each gene, using5 cores...
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=15m 01s
# 287 genes failed the separation test and will be binarized using the fixed cutoff: 0.2
# 3632 genes have lambda2 < 0.03 and will be binarized to 0's.
# Done! Binary assay added to 'assays(sce)$binary'.

og_bin <- sce_p1@assays@data$binary
gs2_bin <- sce2@assays@data$binary

# check the size difference
dim(og_bin)
dim(gs2_bin)

# check how many rows match
str(og_bin)
str(gs2_bin)


# Compare both matrices element-wise to generate a logical matrix.
# R handles the coercion between dense and sparse structures automatically.
mismatch_mat <- og_bin != gs2_bin

# A row matches exactly if it contains zero mismatches across all columns.
exact_row_matches <- rowSums(mismatch_mat) == 0

# Count the total number of identical rows.
total_matching_rows <- sum(exact_row_matches)
print(total_matching_rows)

total_matching_points <- sum(og_bin == gs2_bin)
all_points <- 15048 * 1858
total_matching_points/all_points
# proportion of matching points:
0.4555305

# 1. Flatten both matrices into 1D vectors
# Note: as.matrix() is required to expand the sparse dgCMatrix before vectorising
vec_og <- as.vector(og_bin)
vec_gs2 <- as.vector(as.matrix(gs2_bin))

# 2. Calculate the Phi Coefficient (Pearson Correlation)
# This evaluates the linear dependence between the two vectors
phi_coefficient <- cor(vec_og, vec_gs2)
print(phi_coefficient)

# 3. Calculate Cohen's Kappa
# Generate a contingency table (confusion matrix) to count overlaps
conf_table <- table(Original = vec_og, Refactored = vec_gs2)

# Calculate total elements
total_n <- sum(conf_table)

# Calculate observed agreement (p_o)
p_o <- (conf_table["0", "0"] + conf_table["1", "1"]) / total_n

# Calculate expected chance agreement (p_e)
p_e0 <- (sum(conf_table["0", ]) / total_n) * (sum(conf_table[, "0"]) / total_n)
p_e1 <- (sum(conf_table["1", ]) / total_n) * (sum(conf_table[, "1"]) / total_n)
p_e <- p_e0 + p_e1

# Calculate Kappa
kappa_value <- (p_o - p_e) / (1 - p_e)
print(kappa_value)





###
print(".")
