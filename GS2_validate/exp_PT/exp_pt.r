# LDH1A - experiment:

# I want to know if the incorrect binarisation by the refactored funciton would be better or worse than correct binarisation. 


# assumption, the expression increases throughout the trajectory. (test this by showing the expression across pseudotime )

# if true then if you were to binarise it properly, (mu1 = 0 )
# then you would see the switching point shift earlier in the trajectory.


# maybe this can be done by subjective visualisaiton of exp over pseudotime. 
# if I can identify a clear switching point I can see if at what pseudotime it is. 
# And more importantlu whether it is from 0 to expressed, or from expressed to mroe expressed. 



# library
library(dplyr)
library(tidyr)
library(ggplot2)
library(SingleCellExperiment)


# load the data from the either the original or refactored implementation.
# it doesnt matter because were using expression and pseudotime, which are the same in both implementations.
# sce2
load("/home/mtn1n22/scratch/thesisFigures/GS2_validate/refactored_switching_genes_data.RData")

# use ggplot to visualise the expression of the genes across pseudotime for LDH1A
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/exp_PT/GS2-LDHA-expPT.png", width = 12, height = 6, units = "in", res = 300)
ggplot(data = as.data.frame(t(assay(sce2[rownames(sce2) %in% "LDHA", ], "expdata"))), aes(x = sce2$Pseudotime, y = LDHA)) +
  geom_point(alpha = 0.2) +
  # geom_smooth(method = "loess",span = 2, se = FALSE, color = "blue") +
  labs(title = "Expression of LDHA across Pseudotime",
       x = "Pseudotime",
       y = "Expression") +
  theme_minimal()
dev.off()

# use ggplot to visualise the expression of the genes across pseudotime for PRELID1
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/exp_PT/PRELID1-expPT.png", width = 12, height = 6, units = "in", res = 300)
ggplot(data = as.data.frame(t(assay(sce2[rownames(sce2) %in% "PRELID1", ], "expdata"))), aes(x = sce2$Pseudotime, y = PRELID1)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess",span = 2, se = FALSE, color = "blue") +
  labs(title = "Expression of PRELID1 across Pseudotime",
       x = "Pseudotime",
       y = "Expression") +
  theme_minimal()
dev.off()

# use ggplot to visualise the expression of the genes across pseudotime for TNNI1
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/exp_PT/TNNI1-expPT.png", width = 12, height = 6, units = "in", res = 300)
ggplot(data = as.data.frame(t(assay(sce2[rownames(sce2) %in% "TNNI1", ], "expdata"))), aes(x = sce2$Pseudotime, y = TNNI1)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess",span = 2, se = FALSE, color = "blue") +
  labs(title = "Expression of TNNI1 across Pseudotime",
       x = "Pseudotime",
       y = "Expression") +
  theme_minimal()
dev.off()


# visualise the total expresion of all genes across pseudotime in a single plot
# the data needed needs to be a data frame with columns for pseudotime and total expression 
# therefore we need to calculte the total expression for each cell.

# first calculate total expression for each cell
total_expression <- colSums(assay(sce2, "expdata"))

# create a data frame with pseudotime and total expression
total_expression_df <- data.frame(Pseudotime = sce2$Pseudotime, TotalExpression = total_expression)  

# plot the total expression across pseudotime
png("/home/mtn1n22/scratch/thesisFigures/GS2_validate/exp_PT/TotalExpression-expPT.png", width = 12, height = 6, units = "in", res = 300)
ggplot(data = total_expression_df, aes(x = Pseudotime, y = TotalExpression)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess",span = 2, se = FALSE, color = "", alpha=0.5) +
  labs(title = "Total Expression across Pseudotime",
       x = "Pseudotime",
       y = "Total Expression") +
  theme_minimal()
dev.off()