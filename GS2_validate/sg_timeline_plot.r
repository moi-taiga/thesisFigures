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
png(filename = "tmp.png", width = 10, height = 8, units = "in", res = 300)
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