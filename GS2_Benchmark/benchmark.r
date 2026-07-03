# this script will plot the runtime of the different functions
library(ggplot2)

# load data
runtimes <- read.csv("runtimes.csv")

# runtimes
str(runtimes)
# 'data.frame':   541 obs. of  13 variables:
#  $ Function          : chr  "binarize_exp(current_sce, ncores = 1)" "binarize_exp(current_sce, ncores = 3)" "binarize_exp(current_sce, ncores = 8)" "binarize_exp(current_sce, ncores = 16)" ...
#  $ N_Cores           : int  1 3 8 16 24 32 64 NA NA 1 ...
#  $ Fix_Cutoff        : logi  NA NA NA NA NA NA ...
#  $ Downsample        : logi  NA NA NA NA NA NA ...
#  $ N_Genes           : Factor w/ 8 levels "14759","14826",..: 3 3 3 3 3 3 3 3 3 3 ...
#  $ N_Cells           : Factor w/ 10 levels "1858","3716",..: 1 1 1 1 1 1 1 1 1 1 ...
#  $ User_Time         : logi  NA NA NA NA NA NA ...
#  $ System_Time       : logi  NA NA NA NA NA NA ...
#  $ Elapsed_Time      : num  219.6 87.1 44.4 30.2 26.5 ...
#  $ Total_RAM_Used_MiB: num  215 215 215 215 215 ...
#  $ Peak_RAM_Used_MiB : num  4325 2997 2996 2998 2995 ...
#  $ fun               : Factor w/ 6 levels "glm","glm2","glm_parallel",..: 4 4 4 4 4 4 4 1 1 6 ...
#  $ Peak_RAM_Used_GiB : num  4.22 2.93 2.93 2.93 2.93 ...

# make a column called fun which contains short function name without the parameters, for example "bin_exp" instead of "binarize_exp(current_sce, ncores = 1)"
runtimes$fun <- gsub("\\(.*\\)", "", runtimes$Function)
runtimes$fun <- gsub("binarize_exp_pbmclapply", "refactored", runtimes$fun)
runtimes$fun <- gsub("binarize3_exp", "Arefactored", runtimes$fun)
runtimes$fun <- gsub("binarize_exp", "original", runtimes$fun)

# convert Peak_RAM_Used_MiB to Peak_RAM_Used_GiB
runtimes$Peak_RAM_Used_GiB <- runtimes$Peak_RAM_Used_MiB / 1024

# make N_Genes and N_Cells and fun catagorical variables
runtimes$fun <- as.factor(runtimes$fun)

# plot how fun: original and refactored 's elapsed time changes with number of cores
tmp <- runtimes[runtimes$fun %in% c("original", "refactored"), ]

# only keep where N_Cells is 1858, 9290, 18580
tmp <- tmp[tmp$N_Cells %in% c(1858, 9290, 18580), ]

# set N cores to numeric so x axis is continuous
tmp$N_Cores <- as.numeric(as.character(tmp$N_Cores))
# set ncells to catagorical so shape is discrete
tmp$N_Cells <- as.factor(tmp$N_Cells)

# remove rows where N_Cores is 64
tmp <- tmp[tmp$N_Cores != 64, ]
#remove where n cells is NA
tmp <- tmp[!is.na(tmp$N_Cells), ]
png(filename = "binarize_runtime_VS_cpu.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cores, y = Elapsed_Time, color=fun, shape=N_Cells)) +
    geom_point(size = 4) +
    geom_line() +
    scale_x_continuous(breaks = c(1, 3, 8, 16, 24, 32)) +
    labs(color = "Function", 
         title = "Binarize_exp() - Elapsed Time vs Number of Cores",
             x = "Number of Cores",
             y = "Elapsed Time (seconds)") +
    theme_minimal()
dev.off()



# this time plot RAM*Ncores
# multiplying by N_Cores gives an estimate of the total memory used across all cores,
# this is a pessamistic estimate as mcapply will not coppy the entire object to each core.

# return n cores to numeric
tmp$N_Cores <- as.numeric(as.character(tmp$N_Cores))

tmp$RAM_Ncores <- tmp$Peak_RAM_Used_GiB * tmp$N_Cores
png(filename = "binarize_memory_Ncores_VS_cpu.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cores, y = RAM_Ncores, color=fun, shape = N_Cells)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, aes(group = interaction(fun, N_Cells))) +
  geom_hline(yintercept=240, linetype="dashed", color="grey") + 
  labs(color = "Function",
       title = "Binarize_exp() - Memory Usage vs Number of Cores",
       x = "Number of Cores",
       y = "Peak RAM * N_Cores") +
    theme_minimal()
dev.off()



# plot the elapsed time where fix_cutoff is TRUE
# use N cells as x axis and color bin and bin3
tmp <- runtimes[runtimes$fun %in% c("original", "Arefactored") & runtimes$Fix_Cutoff == TRUE, ]
tmp$fun <- gsub("Arefactored", "refactored", tmp$fun)
# set N-cells to numeric so x axis is continuous
tmp$N_Cells <- as.numeric(as.character(tmp$N_Cells))

# remove rows where cutoff is NA
tmp <- tmp[!is.na(tmp$Fix_Cutoff), ]
# remove the rows where N_Cells is 185800
tmp <- tmp[tmp$N_Cells != 185800, ]

png(filename = "binarize_runtime_VS_Ncells_fix_cutoff.png", width = 10, height = 8, units = "in", res = 300)  
ggplot(tmp, aes(x = N_Cells, y = Elapsed_Time, color=fun)) +
  geom_point(size = 4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 9), se = FALSE, aes(group = fun)) +
  labs(color = "Function",
       title = "Binarize_exp(fix_cutoff = TRUE) - Elapsed Time vs Number of Cells",
       x = "Number of Cells",
       y = "Elapsed Time (seconds)") +
  theme_minimal()
dev.off()

# plot the peak RAM usage where fix_cutoff is TRUE
png(filename = "binarize_memory_VS_Ncells_fix_cutoff.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cells, y = Peak_RAM_Used_GiB, color=fun)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, aes(group = fun)) +
  labs(color = "Function",
       title = "Binarize_exp(fix_cutoff = TRUE) - Peak RAM Used vs Number of Cells",
       x = "Number of Cells",
       y = "Peak RAM Used (GiB)") +
  theme_minimal()
dev.off()


# plot the elapsed time of glm and glm2 

# subset the data to only include glm and glm2
GLMruntimes <- runtimes[runtimes$fun %in% c("find_switch_logistic_fastglm", "find_switch_logistic_fastglm_parallel"), ]

GLMruntimes$fun <- gsub("find_switch_logistic_fastglm_parallel", "refactored", GLMruntimes$fun)
GLMruntimes$fun <- gsub("find_switch_logistic_fastglm", "original", GLMruntimes$fun)


# remove rows where N_Cells is NA
GLMruntimes <- GLMruntimes[!is.na(GLMruntimes$N_Cells), ]
# make sure n cells is numeric
GLMruntimes$N_Cells <- as.numeric(as.character(GLMruntimes$N_Cells))
# remove n cells above 10000
GLMruntimes <- GLMruntimes[GLMruntimes$N_Cells <= 10000, ]

# set n cores to categorical so shape is discrete
GLMruntimes$N_Cores <- as.factor(GLMruntimes$N_Cores)

# where N_Cores is NA, set to 1
GLMruntimes$N_Cores[is.na(GLMruntimes$N_Cores)] <- 1

# filter for downsample TRUE
GLMruntimes <- GLMruntimes[GLMruntimes$Downsample == FALSE, ]
# remove rows where N_Cores is 2
GLMruntimes <- GLMruntimes[GLMruntimes$N_Cores != 2, ]


# x axis is N cells and color is fun
png(filename = "glm_runtime_VS_Ncells.png", width = 10, height = 8, units = "in", res = 300)
ggplot(GLMruntimes, aes(x = N_Cells, y = Elapsed_Time, color=fun, shape= N_Cores)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, aes(group = interaction(fun, N_Cores))) +
  labs(title = "find_switch_logistic_fastglm() - Elapsed Time vs Number of Cells",
       x = "Number of Cells",
       y = "Elapsed Time (seconds)") + 
    theme_minimal()
dev.off()

#set N_Cores to numeric so we can multiply by it
GLMruntimes$N_Cores <- as.numeric(as.character(GLMruntimes$N_Cores))
# make a new column called Total_RAM_ESTIMATE_MiB which is Total_RAM_Used_GiB * N_Cores
GLMruntimes$Total_RAM_ESTIMATE_MiB <- GLMruntimes$Peak_RAM_Used_MiB * GLMruntimes$N_Cores

# set n cores back to factor
GLMruntimes$N_Cores <- as.factor(GLMruntimes$N_Cores)
# plot the memory usage of glm and glm2 
png(filename = "glm_memory_VS_Ncells.png", width = 10, height = 8, units = "in", res = 300)
ggplot(GLMruntimes, aes(x = N_Cells, y = Total_RAM_ESTIMATE_MiB, color=fun, shape= N_Cores)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm",  se = FALSE, aes(group = interaction(fun, N_Cores))) +
  labs(title = "find_switch_logistic_fastglm() - Peak RAM Used vs Number of Cells",
       x = "Number of Cells",
       y = "Peak Total_RAM_ESTIMATE_MiB Used (MiB)") +
    theme_minimal()
dev.off()



# plot the memory usage of glm and glm2 
png(filename = "tmp.png", width = 10, height = 8, units = "in", res = 300)
ggplot(GLMruntimes, aes(x = N_Cells, y = Peak_RAM_Used_MiB, color=fun, shape= N_Cores)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm",  se = FALSE, aes(group = interaction(fun, N_Cores))) +
  labs(title = "find_switch_logistic_fastglm() - Peak RAM Used vs Number of Cells",
       x = "Number of Cells",
       y = "Peak RAM Used (MiB)") +
    theme_minimal()
dev.off()
