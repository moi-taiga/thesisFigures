# this script will plot the runtime of the different functions
library(ggplot2)

# load data
runtimes <- read.csv("runtimes.csv")

# runtimes
str(runtimes)
# 'data.frame':   26 obs. of  12 variables:
#  $ Function          : chr  "binarize_exp(current_sce, ncores = 1)" "binarize_exp(current_sce, ncores = 3)" "binarize_exp(current_sce, ncores = 8)" "binarize_exp(current_sce, ncores = 16)" ...
#  $ N_Cores           : int  1 3 8 16 24 32 64 NA NA NA ...
#  $ Fix_Cutoff        : logi  NA NA NA NA NA NA ...
#  $ Downsample        : logi  NA NA NA NA NA NA ...
#  $ N_Genes           : Factor w/ 2 levels "15048","30096": 1 1 1 1 1 1 1 1 1 1 ...
#  $ N_Cells           : Factor w/ 2 levels "1858","3716": 1 1 1 1 1 1 1 1 1 1 ...
#  $ User_Time         : logi  NA NA NA NA NA NA ...
#  $ System_Time       : logi  NA NA NA NA NA NA ...
#  $ Elapsed_Time      : num  219.6 87.1 44.4 30.2 26.5 ...
#  $ Total_RAM_Used_MiB: num  215 215 215 215 215 ...
#  $ Peak_RAM_Used_MiB : num  4325 2997 2996 2998 2995 ...
#  $ fun               : Factor w/ 4 levels "bin","bin3","glm",..: 1 1 1 1 1 1 1 1 3 3 ...

# make a column called fun which contains short function name without the parameters, for example "bin_exp" instead of "binarize_exp(current_sce, ncores = 1)"
runtimes$fun <- gsub("\\(.*\\)", "", runtimes$Function)
runtimes$fun <- gsub("binarize_exp", "original", runtimes$fun)
runtimes$fun <- gsub("binarize3_exp", "refactored", runtimes$fun)
# convert find_switch_logistic_fastglm to glm
# and find_switch_logistic_fastglm2 to glm2
runtimes$fun <- gsub("find_switch_logistic_fastglm", "glm", runtimes$fun)
runtimes$fun <- gsub("find_switch_logistic_fastglm2", "glm2", runtimes$fun)

# convert Peak_RAM_Used_MiB to Peak_RAM_Used_GiB
runtimes$Peak_RAM_Used_GiB <- runtimes$Peak_RAM_Used_MiB / 1024


# make N_Genes and N_Cells and fun catagorical variables
runtimes$N_Genes <- as.factor(runtimes$N_Genes)
runtimes$N_Cells <- as.factor(runtimes$N_Cells)
runtimes$fun <- as.factor(runtimes$fun)

# plot how fun: original and refactored 's elapsed time changes with number of cores
tmp <- runtimes[runtimes$fun %in% c("original", "refactored"), ]

# tmp <- tmp[tmp$N_Cells == 1858, ]
# remove rows where N_Cores is 64
tmp <- tmp[tmp$N_Cores != 64, ]
#remove where n cells is NA
tmp <- tmp[!is.na(tmp$N_Cells), ]
png(filename = "original_refactored_runtime_VS_cpu.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cores, y = Elapsed_Time, color=fun, shape=N_Cells)) +
    geom_point() +
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

tmp$RAM_Ncores <- tmp$Peak_RAM_Used_GiB * tmp$N_Cores
png(filename = "original_refactored_memory_Ncores_VS_cpu.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cores, y = RAM_Ncores, color=fun, shape = N_Cells)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = c(1, 3, 8, 16, 24, 32, 64)) +
  geom_hline(yintercept=240, linetype="dashed", color="grey") + 
  labs(title = "Memory estimate vs Number of Cores for bin and bin3. \nDashed line at 240 indicates maximum per node",
       x = "Number of Cores",
       y = "memory * N_Cores") +
    theme_minimal()
dev.off()


# plot the elapsed time of bin nd bin3 where fix_cutoff is TRUE
# use N cells as x axis and color bin and bin3
tmp <- runtimes[runtimes$fun %in% c("bin", "bin3") & runtimes$Fix_Cutoff == TRUE, ]
# remove rows where cutoff is NA
tmp <- tmp[!is.na(tmp$Fix_Cutoff), ]
png(filename = "bin_bin3_runtime_VS_Ncells_fix_cutoff.png", width = 10, height = 8, units = "in", res = 300)  
ggplot(tmp, aes(x = N_Cells, y = Elapsed_Time, color=fun)) +
  geom_point() +
  geom_line(aes(group=fun)) +
  labs(title = "Elapsed Time vs Number of Cells for bin and bin3 where fix_cutoff is TRUE.",
       x = "Number of Cells",
       y = "Elapsed Time (seconds)") +
  theme_minimal()
dev.off()


# plot the elapsed time of glm and glm2 
# sahpe by downsample 
# x axis is N cells and color is glm and glm2
tmp <- runtimes[runtimes$fun %in% c("glm", "glm2"), ]
png(filename = "glm_glm2_runtime_VS_Ncells_downsample.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cells, y = Elapsed_Time, color=fun, shape=Downsample)) +
  geom_point(size = 3) +
  geom_line(aes(group=Downsample)) +
  labs(title = "Elapsed Time vs Number of Cells for glm and glm2.",
       x = "Number of Cells",
       y = "Elapsed Time (seconds)") + 
    theme_minimal()
dev.off()

# plot the memory usage of glm and glm2 
tmp <- runtimes[runtimes$fun %in% c("glm", "glm2"), ]
png(filename = "glm_glm2_memory_VS_Ncells_downsample.png", width = 10, height = 8, units = "in", res = 300)
ggplot(tmp, aes(x = N_Cells, y = Peak_RAM_Used_GiB, color=fun, shape=Downsample)) +
  geom_point(size = 3) +
  geom_line(aes(group=Downsample)) +
  labs(title = "Peak RAM Used vs Number of Cells for glm and glm2.",
       x = "Number of Cells",
       y = "Peak RAM Used (GiB)") +
    theme_minimal()
dev.off()