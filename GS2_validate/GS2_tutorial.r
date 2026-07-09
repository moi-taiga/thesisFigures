## Input log-normalized gene expression, Monocle2 pseudo-time and dimensionality reduction
## Path1 containing cells in states 3,2,1
sce2 <- convert_monocle2(monocle2_obj = cardiac_monocle2,
                           states = c(3,2,1), expdata = logexpdata)

### binarize gene expression using gene-specific thresholds
sce2 <- gs2_bin(sce2, ncores = 7)

## fit logistic regression and find the switching pseudo-time point for each gene
## with downsampling. This step takes less than 1 mins
# library fas
library(fastglm)
sce2 <- gs2_glm(sce2, downsample = TRUE, ncores = 3)
sg2_allgenes <- filter_switchgenes(sce2, allgenes = TRUE, r2cutoff = 0.00)


# save the timeline plot of all genes
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/refactored_switching_genes_timeline.png", width = 10, height = 8, units = "in", res = 300)
plot_timeline_ggplot(sg2_allgenes, timedata = sce2$Pseudotime, txtsize = 3)
dev.off()

# save the sce2 and sg2object
save(sce2, sg2_allgenes, file = "/home/mtn1n22/scratch/thesisFigures/GS2_validate/refactored_switching_genes_data.RData")
