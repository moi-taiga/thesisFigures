##### old bin ###########

old_bin <- function(sce, fix_cutoff = FALSE, binarize_cutoff = 0.2, ncores = 3) {
  # calculate zero percentage
  zerop_g <- c()
  expdata <- assays(sce)$expdata
  for (i in 1:nrow(expdata)) {
    zp <- length(which(expdata[i, ] == 0))/ncol(expdata)
    zerop_g <- c(zerop_g, zp)
  }

  if (fix_cutoff == TRUE) {
    expdata <- assays(sce)$expdata
    is.na(expdata) <- assays(sce)$expdata == 0
    exp_reduced_binary <- as.matrix((expdata > binarize_cutoff) + 0)
    exp_reduced_binary[is.na(exp_reduced_binary)] = 0
    assays(sce)$binary <- exp_reduced_binary
    oupBinary <- data.frame(geneID = rownames(sce),
                            zerop_gene = zerop_g,
                            passBinary = TRUE)
    rowData(sce) <- oupBinary
  } else {
    expdata <- assays(sce)$expdata
    # Add gaussian noise to gene expression matrix
    # Here we use a sd of 0.1
    LogCountsadd = expdata + matrix(rnorm(nrow(expdata)*ncol(expdata),
                                          mean = 0, sd = 0.1),
                                    nrow(expdata), ncol(expdata))
    # Start fitting mixture models for each gene
    oupBinary = do.call(
      rbind, mclapply(rownames(LogCountsadd), function(iGene){
        set.seed(42)   # Set seed for consistency
        tmpMix = normalmixEM(LogCountsadd[iGene, ], k = 2)
        if (tmpMix$mu[1] < tmpMix$mu[2]) {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[1],
                              mu2 = tmpMix$mu[2],
                              sigma1 = tmpMix$sigma[1],
                              sigma2 = tmpMix$sigma[2],
                              lambda1 = tmpMix$lambda[1],
                              lambda2 = tmpMix$lambda[2],
                              loglik = tmpMix$loglik)
        } else {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[2],
                              mu2 = tmpMix$mu[1],
                              sigma1 = tmpMix$sigma[2],
                              sigma2 = tmpMix$sigma[1],
                              lambda1 = tmpMix$lambda[2],
                              lambda2 = tmpMix$lambda[1],
                              loglik = tmpMix$loglik)
        }
        return(tmpOup)
      }, mc.cores = ncores))

    # Check if non-bimodal genes
    oupBinary$passBinary = TRUE
    oupBinary[oupBinary$lambda1 < 0.1, ]$passBinary = FALSE
    oupBinary[oupBinary$lambda2 < 0.1, ]$passBinary = FALSE
    oupBinary[(oupBinary$mu2 - oupBinary$mu1) < (oupBinary$sigma1 + oupBinary$sigma2), ]$passBinary = FALSE
    # table(oupBinary$passBinary)

    # Solve for intersection for remaining genes
    oupBinary$root = -1
    for(iGene in oupBinary[oupBinary$passBinary == TRUE, ]$geneID){
      tmpMix = oupBinary[oupBinary$geneID == iGene, ]
      tmpInt = uniroot(function(x, l1, l2, mu1, mu2, sd1, sd2) {
        dnorm(x, m = mu1, sd = sd1) * l1 -
          dnorm(x, m = mu2, sd = sd2) * l2},
        interval = c(tmpMix$mu1,tmpMix$mu2),
        l1 = tmpMix$lambda1, mu1 = tmpMix$mu1, sd1 = tmpMix$sigma1,
        l2 = tmpMix$lambda2, mu2 = tmpMix$mu2, sd2 = tmpMix$sigma2)
      oupBinary[oupBinary$geneID == iGene, ]$root = tmpInt$root
    }
    # Binarize expression
    binLogCounts = expdata[oupBinary$geneID,]
    binLogCounts = t(scale(t(binLogCounts), scale = FALSE,
                           center = oupBinary$root))
    binLogCounts[binLogCounts >= 0] = 1
    binLogCounts[binLogCounts < 0] = 0
    assays(sce)$binary <- binLogCounts

    oupBinary$zerop_gene <- zerop_g
    rowData(sce) <- oupBinary
  }
  return(sce)
}


########## old glm ############


downsample_zeros <- function(glmdata, ratio_ds = 0.7) {
  p = as.numeric(ratio_ds)
  set.seed(42)   # Set seed for consistency
  downsample <- sample(which(glmdata$expvalue == 0), length(which(glmdata$expvalue == 0)) - round(sum(glmdata$expvalue != 0) * p/(1 - p)))
  if (length(downsample) > 0) {
    subdata <- glmdata[-downsample, ]
  } else {subdata <- glmdata}
  return(subdata)
}


old_glm <- function(sce, downsample = FALSE, ds_cutoff = 0.7, zero_ratio = 0.7,
                                         sig_FDR = 0.05, show_warnings = TRUE) {
  binarydata <- assays(sce)$binary
  expdata <- assays(sce)$expdata
  binarydata <- binarydata[which(rowData(sce)$passBinary == TRUE), ]
  expdata <- expdata[which(rowData(sce)$passBinary == TRUE), ]
  genes <- rowData(sce)[which(rowData(sce)$passBinary == TRUE), ]
  timedata <- sce$Pseudotime
  pvalues <- binarydata[, 1]
  pseudoR2s <- binarydata[, 1]
  estimates <- binarydata[, 1]
  switch_at_time <- binarydata[, 1]
  prd_quality <- binarydata[, 1]
  CI <- binarydata[, 1]

  for (i in 1:nrow(binarydata)) {
    glmdata <- cbind(State = as.numeric(binarydata[i, ]), expvalue = as.numeric(expdata[i, ]),
                     timedata = sce$Pseudotime)
    glmdata <- as.data.frame(glmdata)

    if (downsample == TRUE & round(genes$zerop_gene[i],3) > ds_cutoff) {
      glmdata <- downsample_zeros(glmdata, ratio_ds = zero_ratio)
    }

    if (show_warnings == TRUE) {
      glm_results <- fastglm(x = model.matrix(State ~ timedata, data = glmdata),
                             y = glmdata$State, family = binomial(link = "logit"))
    } else {
      glm_results <-suppressWarnings(fastglm(x = model.matrix(State ~ timedata, data = glmdata),
                                             y = glmdata$State, family = binomial(link = "logit")))
    }
    pvalues[i] <- coef(summary(glm_results))[, 4][2]
    ll.null <- glm_results$null.deviance/-2
    ll.proposed <- glm_results$deviance/-2
    # McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
    pseudoR2s[i] <- (ll.null - ll.proposed)/ll.null
    estimates[i] <- coef(summary(glm_results))[, 1][2]
    # p=0.5
    switch_at_time[i] <- (log(0.5/(1 - 0.5)) - coef(glm_results)[1])/coef(glm_results)[2]
    if (switch_at_time[i] >= max(glmdata$timedata)) {
      switch_at_time[i] = max(glmdata$timedata)
      prd_quality[i] = 0
    } else {
      prd_quality[i] = 1
    }
    if (switch_at_time[i] <= min(glmdata$timedata)) {
      switch_at_time[i] = min(glmdata$timedata)
      prd_quality[i] = 0
    }
    se <- summary(glm_results)$coefficients[, 2]
    CI[i] <- sqrt((se[1]*1.96/coef(glm_results)[1])^2 + (se[2]*1.96/coef(glm_results)[2])^2)*
      abs(coef(glm_results)[1]/coef(glm_results)[2])
    remove(glm_results)
  }

  result_switch <- cbind(switch_at_time, CI, pvalues, pseudoR2s, estimates, prd_quality)
  rownames(result_switch) <- rownames(binarydata)
  result_switch <- as.data.frame(result_switch)
  result_switch$direction <- "up"
  result_switch[result_switch$estimates < 0, ]$direction <- "down"
  result_switch$FDR <- p.adjust(result_switch$pvalues, method = "BH")
  steptime <- (max(timedata) - min(timedata))/100
  result_switch$switch_at_timeidx <- round((result_switch$switch_at_time - min(timedata))/steptime)

  # process_resultswitch --------------------------------------------------------------------------
  # check significance FDR < sig_FDR
  if (max(result_switch$FDR) > sig_FDR) {
    result_switch[result_switch$FDR > sig_FDR, ]$prd_quality <- 0
  }

  geneinfo <- merge(rowData(sce), result_switch, by=0, all=TRUE)[,-1] #[,1:11]
  rownames(geneinfo) <- geneinfo$geneID
  # all(rownames(geneinfo) == rownames(sce))
  geneinfo <- geneinfo[rownames(sce), ]
  rowData(sce) <- geneinfo
  return(sce)
}