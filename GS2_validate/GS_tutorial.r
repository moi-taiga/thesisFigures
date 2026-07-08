## load libraries
library(GeneSwitches)
library(SingleCellExperiment)

## Download example files to current directory
#get_example_inputData()
## Load input data log-normalized gene expression
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/logexpdata.RData")
## Load Monocle2 object with pseudo-time and dimensionality reduction
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/cardiac_monocle2.RData")

## Input log-normalized gene expression, Monocle2 pseudo-time and dimensionality reduction
## Path1 containing cells in states 3,2,1
sce_p1 <- convert_monocle2(monocle2_obj = cardiac_monocle2,
                           states = c(3,2,1), expdata = logexpdata)

### binarize gene expression using gene-specific thresholds
library(mclapply)
library(parallel)
library(mixtools)
sce_p1 <- old_bin(sce_p1, ncores = 4)

## fit logistic regression and find the switching pseudo-time point for each gene
## with downsampling. This step takes less than 1 mins
library(fastglm)
sce_p1 <- old_glm(sce_p1, downsample = TRUE, show_warning = FALSE)
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.00)

#save the timeline plot of all genes
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/old_switching_genes_timeline.png", width = 10, height = 8, units = "in", res = 300)
plot_timeline_ggplot(sg_allgenes, timedata = sce_p1$Pseudotime, txtsize = 3)
dev.off()

# save the sce_p1 and sgobject
save(sce_p1, sg_allgenes, file = "/home/mtn1n22/scratch/thesisFigures/GS2_validate/old_switching_genes_data.RData")
