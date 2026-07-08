# new gs2 funcitons

#### binarize ####
#' @title Binarize gene expression
#'
#' @description This function generates on/off binarized data for gene expression
#'
#' @param sce SingleCellExperiment
#' @param fix_cutoff Logical. if use fixed global cutoff for binarization, default FALSE
#' @param binarize_cutoff fixed global cutoff for binarization, default 0.2
#' @param ncores number of cores
#' @param gaussian_weight_cutoff float if the second gaussian has a weight (lambda) less than this cutoff assume the gene is not expressed.
#' @return A SingleCellExperiment object with an added binary assay in `assays(sce)$binary` and updated `rowData(sce)` containing binarization metadata.
#'
#' @import parallel
#' @importFrom pbapply pblapply
#' @importFrom mixtools normalmixEM
#' @importFrom Matrix rowSums colSums t Matrix
#' @export
#'
gs2_bin <- function(sce, fix_cutoff = FALSE, binarize_cutoff = 0.2, ncores = 3, gaussian_weight_cutoff = 0.03) {

  # Check if binarize_cutoff is valid -must be greater than 0 and numeric
  if (!is.numeric(binarize_cutoff) || binarize_cutoff < 0) {
    stop("The 'binarize_cutoff' must be a numeric value greater than or equal to 0.")
  }

  # Check if binarized assay already exists to avoid overwriting
  if ("binary" %in% names(assays(sce))) {
    stop("The 'binary' assay already exists in the SingleCellExperiment object. \nPlease remove or rename it before running this function.")
  }

  # Check if expdata exists
  if (is.null(assays(sce)$expdata)) {
    stop("Expression data 'expdata' not found in the SingleCellExperiment object assays.")
  }

  # load in the expression data from the sce object
  expdata <- assays(sce)$expdata

  # ensure expdata is a sparse matrix for downstream efficiency
  if (!inherits(expdata, "dgCMatrix")) {
    message("assays(sce)$expdata is not sparse (dgCMatrix).\nCoercing expression data to dgCMatrix for efficiency...")
    expdata <- as(expdata, "dgCMatrix")
  }

  # calculate the percentage of zeros for each gene
  # this variable is not used in the current binarization funciton but used downstream?
  # assumes expdata is a sparse matrix
  # uses Matrix package to calculate the percentage of zeros for each gene
  zerop_g <- 1 - (Matrix::rowSums(expdata != 0) / ncol(expdata))

  if (fix_cutoff == TRUE) {
    exp_reduced_binary <- (expdata > binarize_cutoff) + 0 # Remains sparse 'dgCMatrix'
    assays(sce)$binary <- exp_reduced_binary
    oupBinary <- data.frame(geneID = rownames(sce),
                            zerop_gene = zerop_g,
                            passBinary = FALSE) # as is has not passed any bimodality tests.
    rowData(sce) <- oupBinary
  } else {
    # transpose the expdata to have genes as columns and cells as rows for easier processing in the parallel loop
    expdata_t <- Matrix::t(expdata)

    # before fitting the mixture models, filter out genes that are expressed in less than 5 cells,
    # as these are likely to fail the bimodality tests and cause long runtimes due to non-convergence of the normalmixEM algorithm.
    # these will later be binarized to all 0's, which is likely more accurate than trying to fit a bimodal distribution to a gene that is only expressed in a few cells.
    n_cells_expressed <- Matrix::colSums(expdata_t > 0)
    low5_exp_genes <- names(n_cells_expressed)[n_cells_expressed <= 5]
    if (length(low5_exp_genes) > 0) {
      message(paste(length(low5_exp_genes), "genes are expressed in 5 or fewer cells and will be binarized to 0's."))
    }

    # filter out these genes from expdata_t for the mixture model fitting
    expdata_t <- expdata_t[, !(colnames(expdata_t) %in% low5_exp_genes)]

    # --- MIXTURE MODEL FITTING ---
    message("Fitting mixture models for each gene, using", ncores, " cores...")

    # Define the function to fit mixture models for a single gene
    fit_gene_fun <- function(iGene) {
      # START TRYCATCH to avoid errors from the normalmixEM function which can fail to converge for some genes.
      tryCatch({

        # Extract sparse vector for gene
        raw_counts <- as.numeric(expdata_t[, iGene])

        # Add gaussian noise to the gene expression to smooth the distribution for better fitting of mixture models.
        # # Here we use a sd of 0.1
        set.seed(42)   # Set seed for consistency
        noisy_counts <- raw_counts + rnorm(length(raw_counts), mean = 0, sd = 0.1)


        # fit mixture model with two components per gene
        # stricter max itterations and restarts to avoid long runtimes on genes that fail to converge.
        # similar logic for the higher epsilon value.
        # silence messages as they clutter the terminal.
        invisible(capture.output({
          tmpMix = mixtools::normalmixEM(noisy_counts,
                                         k = 2,
                                         maxit = 400,
                                         maxrestarts = 15,
                                         epsilon = 1e-6,
                                         verb = FALSE)
        }))


        if (tmpMix$mu[1] < tmpMix$mu[2]) {
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[1],
                              mu2 = tmpMix$mu[2],
                              sigma1 = tmpMix$sigma[1],
                              sigma2 = tmpMix$sigma[2],
                              lambda1 = tmpMix$lambda[1],
                              lambda2 = tmpMix$lambda[2],
                              loglik = tmpMix$loglik,
                              passBinary = TRUE)
        } else { # correct the order of the two gausians.
          tmpOup = data.frame(geneID = iGene,
                              mu1 = tmpMix$mu[2],
                              mu2 = tmpMix$mu[1],
                              sigma1 = tmpMix$sigma[2],
                              sigma2 = tmpMix$sigma[1],
                              lambda1 = tmpMix$lambda[2],
                              lambda2 = tmpMix$lambda[1],
                              loglik = tmpMix$loglik,
                              passBinary = TRUE)
        }
        return(tmpOup)
      }, error = function(e) {
        # ERROR HANDLER: If normalmixEM fails, return NA values for this gene
        message("Error fitting gene: ", iGene)
        # likley due to non-convergence
        # Return a row with NA stats and passBinary = FALSE
        data.frame(geneID = iGene,
                   mu1 = NA, mu2 = NA,
                   sigma1 = NA, sigma2 = NA,
                   lambda1 = NA, lambda2 = NA,
                   loglik = NA,
                   passBinary = FALSE)
      })
      # END TRYCATCH
    }

    # Start fitting mixture models for each gene
    # If on Windows, use pbapply. If on Linux/Mac, use pbmcapply with (more memory efficient).
    if (.Platform$OS.type == "windows") {
      oupBinary_list <- pbapply::pblapply(colnames(expdata_t), fit_gene_fun, cl = ncores)
    } else {
      oupBinary_list <- pbmcapply::pbmclapply(colnames(expdata_t), fit_gene_fun, mc.cores = ncores)
    }

    oupBinary <- do.call(rbind, oupBinary_list)

    # report number of genes which normalmixEM failed to fit a model
    failed_fit_count <- sum(oupBinary$passBinary == FALSE)
    if (failed_fit_count > 0) {
      message(paste(failed_fit_count, "genes failed to fit a mixture model and will be binarized using the fixed cutoff:", binarize_cutoff))
    }


    # --- Biomodality Tests ---

    # Check for genes that do not pass the bimodality tests

    # check for overlap between the two gausians,
    # if the distance between the means is less than the sum of the SDs, mark as non-bimodal
    # assume these are continuous genes that do not have a clear switch point and are not good candidates for binarization.

    # using which(oupBinary$passBinary) to ensure we don't calculate on NAs
    pass_idx <- which(oupBinary$passBinary) # Only look at currently passing genes
    diff_mu <- oupBinary$mu2[pass_idx] - oupBinary$mu1[pass_idx]
    sum_sigma <- oupBinary$sigma1[pass_idx] + oupBinary$sigma2[pass_idx]
    fail_sep <- pass_idx[diff_mu < sum_sigma]
    oupBinary$passBinary[fail_sep] <- FALSE
    # table(oupBinary$passBinary)
    # print the number of genes that fail the separation test
    if (length(fail_sep) > 0) {
      message(paste(length(fail_sep), "genes failed the separation test and will be binarized using the fixed cutoff:", binarize_cutoff))
    }

    # We must use 'which' to avoid crashing on the genes that failed above (they introduce NA values)
    # Check if one gaussian dominates the other, if one gaussian has less than x% weight mark as non-bimodal
    # changing this to a variable, default should be more like ~3% based on my testing. - MTN 05/02/26

    # report the number of genes that fail the gaussian weight cutoff for each gaussian.
    idx_l2_too_small <- which(oupBinary$passBinary & oupBinary$lambda2 < gaussian_weight_cutoff)
    if (length(idx_l2_too_small) > 0) {
      message(paste(length(idx_l2_too_small), "genes have lambda2 <", gaussian_weight_cutoff, "and will be binarized to 0's."))
    }
    # Note: the logic above relies on the testing of distance between the means being done first.

    # filter for l1 being too small removed as it would remove housekeeping genes,
    # I cant find many/any scenarios where a gene with a small l1 would be a problem for binarization

    # These genes should be marked as non-bimodal and set to 0 expression in the binarization
    # initialize a new column in oupBinary to store the roots,
    # doing this early alows me to add Inf values for the genes that fail the lambda2 bimodality test.
    oupBinary$root <- NA
    # Inf values mean these genes will be binarized to all 0's
    # only do this for genes where oupBinary$passBinary is TRUE to avoid adding Inf values to genes that already failed the separation test and have passBinary = FALSE.
    # Use the integer indices (idx_l2_too_small) directly; they already filtered for passBinary == TRUE
    oupBinary$root[idx_l2_too_small] <- Inf

    # --- CALCULATE Roots/intersections ---

    # Identify the intersection between the two gausians
    # aka the switching point of a switching gene

    # Identify the indicies of the genes which passed the tests for bimodality
    pass_indices = which(oupBinary$passBinary == TRUE & is.na(oupBinary$root))

    # Itterate through each gene and find the root/intersection
    for(i in pass_indices){
      # select the values for the current gene
      mu1 = oupBinary$mu1[i]
      mu2 = oupBinary$mu2[i]
      sd1 = oupBinary$sigma1[i]
      sd2 = oupBinary$sigma2[i]
      l1 = oupBinary$lambda1[i]
      l2 = oupBinary$lambda2[i]

      # Added tryCatch here because in rare cases the intersection is not between the two means, which causes uniroot to fail.
      tryCatch({
        # calculate the root/intersection between the two gaussians using uniroot
        tmpInt = uniroot(function(x) {
          dnorm(x, mean = mu1, sd = sd1) * l1 -
            dnorm(x, mean = mu2, sd = sd2) * l2},
          interval = c(mu1, mu2))
        # store the root in the dataframe
        oupBinary$root[i] = tmpInt$root

      }, error = function(e) {
        # If uniroot fails, we mark the gene as failed
        oupBinary$passBinary[i] <<- FALSE # Use double arrow to update outside function
      })

    }


    # --- BINARIZATION ---

    # match the order of genes in expdata_t
    # (because oupBinary might be in a different order after the parallel loop)
    rownames(oupBinary) <- oupBinary$geneID
    oupBinary <- oupBinary[colnames(expdata_t), ]

    # Identify genes that failed constraints or model fitting (where root is NA or passBinary is FALSE)
    # Instead of removing them, we assign the fixed binarize_cutoff.
    genes_failed <- which(is.na(oupBinary$root) | oupBinary$passBinary == FALSE)

    if (length(genes_failed) > 0) {
      # Set the root to the manual cutoff for these genes
      oupBinary$root[genes_failed] <- binarize_cutoff
    }


    # Perform Sparse Binarization efficiently
    # binarize based on the root value for each gene,
    # if the expression is greater than the root, it is 1, otherwise 0.

    # Transpose back to Genes x Cells immediately.
    # We do this because we need to compare each row (gene) to a specific threshold.
    mat_genes <- Matrix::t(expdata_t)

    # Force to dgCMatrix to ensure @i and @x slots are accessible
    # (t() can sometimes return dgTMatrix or other types)
    if (!inherits(mat_genes, "dgCMatrix")) {
      mat_genes <- as(mat_genes, "dgCMatrix")
    }

    # Extract thresholds aligned with the rows of mat_genes
    gene_roots <- oupBinary$root

    # We leverage the internal structure of dgCMatrix:
    # @x contains non-zero values, @i contains their 0-based row indices.
    # We only need to update the non-zero values.
    # (Assumption: gene_roots >= 0. Therefore 0 > root is always FALSE, so zeros stay zeros).
    # aka they are already binarized to 0's, so we only need to update the non-zero values to 1's or 0's based on the root comparison.

    # Get 1-based row indices for every non-zero element
    nonZero_row_indices <- mat_genes@i + 1

    # Compare each non-zero value to its specific gene threshold
    # If val > root then expression is 1. If val <= root then expression is 0.
    # Note: If root is Inf, val > Inf is FALSE then expression is 0.

    # Update the matrix values
    mat_genes@x  <- (mat_genes@x > gene_roots[nonZero_row_indices]) + 0

    # Drop explicit zeros (values that were <= root became 0) to keep matrix sparse
    binLogCounts <- Matrix::drop0(mat_genes)

    # add back the genes that were filtered out for being expressed in less than 5 cells, and binarize them to all 0's.
    if (length(low5_exp_genes) > 0) {
      # Create sparse matrix of 0s:
      # with the same number of columns as low5_exp_genes and the same number of rows as expdata_t.
      low5_binary <- Matrix::Matrix(0, nrow = length(low5_exp_genes), ncol = ncol(binLogCounts),
                                    dimnames = list(low5_exp_genes, colnames(binLogCounts)))

      # also add them to the metadata dataframe oupBinary
      # add the skipped genes back to oupBinary
      # with passBinary = TURE (as they will be binarized to all 0's, and should be used for downstream analysis)
      # set the other parameters asside from geneID and passBinary to NA as they were not tested.
      oupSkipped_genes <- data.frame(geneID = low5_exp_genes,
                                     mu1 = NA, mu2 = NA,
                                     sigma1 = NA, sigma2 = NA,
                                     lambda1 = NA, lambda2 = NA,
                                     loglik = NA,
                                     passBinary = FALSE, # As they were not tested for bimodality
                                     root = Inf) # Infinite root means always 0
      rownames(oupSkipped_genes) <- oupSkipped_genes$geneID
      oupBinary <- rbind(oupBinary, oupSkipped_genes)

      # Bind this with the binarized matrix
      # (Genes as Rows as we have transposed back)
      binLogCounts <- rbind(binLogCounts, low5_binary)
    }

    # ensure the order of genes in binLogCounts matches the order of genes in sce
    binLogCounts <- binLogCounts[rownames(sce), ]
    # Store in SingleCellExperiment
    assays(sce)$binary <- binLogCounts

    # Add metdata
    # ensure the order of genes in oupBinary matches the order of genes in sce
    oupBinary <- oupBinary[rownames(sce), ]
    # zerop_g to oupBinary
    oupBinary$zerop_gene <- zerop_g[oupBinary$geneID]
    # oupBinary to rowData of sce
    rowData(sce) <- oupBinary


    message("Done! Binary assay added to 'assays(sce)$binary'.")
    #add a warning if binary includes NA's # They should not include NA's as we have removed the genes that fail the bimodality tests, but this is just a safety check.
    if (any(is.na(assays(sce)$binary))) {
      warning("The binary assay contains NA values, which may cause issues in downstream analysis if not removed.")
    }
  }
  return(sce)
}
#### glm ####
#' @title Fit fast logistic regression and find switching timepoint
#'
#' @description This function fits fast logistic regression and find switching timepoint for each gene
#'
#' @param sce SingleCellExperiment
#' @param downsample Logical. if do random downsampling of zeros
#' @param ds_cutoff only do downsampling if zero percentage is over this cutoff
#' @param zero_ratio downsampling zeros to this proportion
#' @param sig_FDR FDR cut off for significant genes
#' @param ncores Numeric. Number of cores to use for parallel processing. Default is 3.
#' @param EPV Numeric. Minimum number of Events Per Variable (EPV) required for a gene to be included in the regression analysis. Default is 10.
#' @return A SingleCellExperiment object with the regression results added to the rowData.
#'
#' @import fastglm
#' @export
#'
gs2_glm <- function(sce,
                                         downsample = TRUE,
                                         ds_cutoff = 0.7,
                                         zero_ratio = ds_cutoff,
                                         sig_FDR = 0.05,
                                         ncores = 3,
                                         EPV = 10) {
  # CHECK: Ensure the necessary assays are present in the SingleCellExperiment object
  if (!("binary" %in% names(assays(sce))) || !("expdata" %in% names(assays(sce)))) {
    stop("The SingleCellExperiment object must contain 'binary' and 'expdata' assays. Please run 'binarize_exp()' first to create these assays.")
  }

  # load the data binarized by binarize_exp(), and the original expression data
  binarydata <- assays(sce)$binary
  expdata <- assays(sce)$expdata

  # this function assumes the matrices to be sparse
  # ensure expdata is a sparse matrix for downstream efficiency
  if (!inherits(expdata, "dgCMatrix")) {
    message("assays(sce)$expdata is not sparse (dgCMatrix).",
            "Coercing expression data to dgCMatrix for efficiency...")
    expdata <- as(expdata, "dgCMatrix")
  }
  # ensure binarydata is a sparse matrix for downstream efficiency
  if (!inherits(binarydata, "dgCMatrix")) {
    message("assays(sce)$binary is not sparse (dgCMatrix).",
            "Coercing expression data to dgCMatrix for efficiency...")
    binarydata <- as(binarydata, "dgCMatrix")
  }

  # --- TRANSPOSE ---
  # looping through a sparse matrix's columns is faster than looping through rows.
  binary_t <- Matrix::t(binarydata)
  exp_t <- Matrix::t(expdata)

  # --- Load/PRE-CALCULATE CONSTANTS ---
  # initialize vectors to store the results of the regression for each gene
  # these vectors have the same length as the number of genes
  # intialize as NA, for the genes which fail regression.
  # and contan the gene names as rownames
  n_genes <- ncol(binary_t)
  n_cells <- nrow(binary_t)

  # Pre-calculate Time Bounds
  # extract pseudotime data from the sce object
  timedata <- sce$Pseudotime
  # The loop calls min() and max() on 'timedata' repeatedly.
  # Since pseudotime is the same for all genes, calculate this once.
  t_min <- min(timedata)
  t_max <- max(timedata)

  # Pre-calculate the Design Matrix
  # Instead of building this matrix ~20,000 times inside the loop, build it once.
  # This matrix contains the Intercept (col 1) and Pseudotime (col 2).
  design_mat_full <- model.matrix(~ timedata)

  # Pre-calculate Downsampling Trigger
  # Instead of doing the 'round()' and boolean checks inside the loop,
  # create a simple TRUE/FALSE vector for all genes at once.
  if (downsample) {
    should_downsample <- round(rowData(sce)$zerop_gene, 3) > ds_cutoff  # not sure why we need to round here, but keeping consistent with original code -MTN 08/02/26
  } else {
    should_downsample <- logical(n_genes) # Initialize a vector of FALSE values for all genes
  }

  # SPARSE CHECK - Skip genes where the regression is unlikely to work due to too few 1s or too few 0s.
  # Calculate the number of cells expressing each gene (Sum of the column in binary_t)
  gene_counts <- Matrix::colSums(binary_t)
  # Define the threshold
  min_cells_required <- EPV
  # Identify genes to skip:
  # 1. Expressed in < 10 cells (Too few 1s)
  # 2. Expressed in > (Total - 10) cells (Too few 0s)
  skip_gene <- (gene_counts < min_cells_required) | (gene_counts > (n_cells - min_cells_required))
  # report number of skipped genes
  message(sum(skip_gene), " genes will be skipped due to insufficient expression.")
  # "If a gene is expressed in fewer than 10 cells regression results may be unreliable"

  message("Fitting logistic regression models using ", ncores, " cores...")

  # Define the worker function for parallelization
  fit_single_gene <- function(i) {
    # Default return for skipped/failed genes
    fail_result <- list(
      pvalue = NA_real_,
      pseudoR2 = NA_real_,
      estimate = NA_real_,
      switch_at_time = NA_real_,
      prd_quality = 0,
      CI = NA_real_
    )

    # Skip genes that are too sparse for reliable regression
    if (skip_gene[i]) {
      return(fail_result)
    }

    # extract the binary state for the current gene
    gene_exp_vec <- as.numeric(binary_t[, i])

    #
    glm_response_var <- gene_exp_vec
    glm_predictor_var <- design_mat_full

    # --- DOWNSAMPLING ZEROS---
    # downsample zeros if the percentage of zeros is above the cutoff and downsampling is enabled
    # Downsample inside the loop, but only for genes that meet the downsampling criteria.
    # This methods avoids building a new dataframe for every gene, and instead just slices the existing vectors/matrices.
    if (should_downsample[i]) {
      # We only fetch the raw expression data IF we actually need to downsample
      exp_vec <- as.numeric(exp_t[, i])

      # Identify zeros
      zero_indices <- which(exp_vec == 0)
      n_non_zero   <- length(exp_vec) - length(zero_indices)

      # Calculate number of zeros to keep based on the desired ratio
      n_keep <- round(n_non_zero * (zero_ratio / (1 - zero_ratio)))

      # Sample if necessary
      set.seed(42) #
      keep_zeros <- sample(zero_indices, n_keep)

      # Combine indices: (All non-zeros) + (Selected zeros)
      keep_idx <- c(which(exp_vec != 0), keep_zeros)

      # Update our working variables by slicing
      glm_response_var <- gene_exp_vec[keep_idx]
      glm_predictor_var <- design_mat_full[keep_idx, , drop = FALSE]
    }

    # fit a logistic regression model using fastglm,
    # with the binary expression as the response variable and pseudotime as the predictor variable
    # Fit the model using the prepared matrices
    # use trycatch here if errors occur due to non-convergence or other issues with the regression.
    tryCatch({
      fit <- fastglm(x = glm_predictor_var, y = glm_response_var, family = binomial(link = "logit"))
    }, error = function(e) {
      warning(paste("Model fitting failed for gene index", i, ":", e$message))
      return(fail_result)
    })


    # If fit failed (NULL), skip to next gene
    if (is.null(fit)) {
      return(fail_result)
    }

    # --- EXTRACT REGRESSION RESULTS ---
    ## # extract information from the regression results.
    # # calculate the time at which the predicted probability of expression is 0.5 using the coefficients from the logistic regression model.
    # # check if the switching time is within the range of the pseudotime data, and adjust it if it is outside the range.
    # # this could be optional, as in some (rare) downstream applications, it might be useful to know  which genes are predicted to switch before or after the observed pseudotime range.
    # calculate the confidence interval for the switching time using the standard errors of the coefficients from the logistic regression model.

    # Extract Summary ONCE
    sum_fit <- summary(fit)
    coefs   <- coef(sum_fit)

    # Check for valid output (needs at least 2 rows: Intercept + Time)
    if (nrow(coefs) < 2) {
      return(fail_result)
    }

    #Store Basic Stats
    p_val   <- coefs[2, 4] # P-value (Row 2, Col 4)
    est_val <- coefs[2, 1] # Slope (Row 2, Col 1)

    #Pseudo R2
    ll_null     <- fit$null.deviance / -2
    ll_proposed <- fit$deviance / -2
    pseudoR2_val <- (ll_null - ll_proposed) / ll_null

    #Switch Time Calculation
    # Formula: Time = (TargetLogit - Intercept) / Slope
    # Since TargetLogit for 0.5 is 0, Time = -Intercept / Slope
    intercept_val <- coefs[1, 1]
    slope_val     <- coefs[2, 1]

    calc_time <- -intercept_val / slope_val

    switch_time_val <- NA_real_
    quality_val <- 0

    # Quality Control (Using pre-calculated t_min/t_max)
    # if the calculated switching time is outside the observed pseudotime range,
    #   we set the switching time to the nearest bound (t_min or t_max) and mark the prediction quality as 0 (low).
    if (calc_time >= t_max) {
      switch_time_val <- t_max; quality_val <- 0
    } else if (calc_time <= t_min) {
      switch_time_val <- t_min; quality_val <- 0
    } else {
      switch_time_val <- calc_time; quality_val <- 1
    }

    # Confidence Interval (Delta Method)
    # combines the standard errors of the intercept and slope to estimate the uncertainty in the switching time.
    # if statement to check if calc_time is NA (which can happen if slope_val is 0)
    # or if the intercept is very close to zero (which can lead to errors in the CI calculation).
    # If either of these conditions is true, we set the CI to NA to avoid misleading results.
    CI_val <- NA_real_
    if (!is.na(calc_time) && abs(intercept_val) > 1e-10) {
      se_intercept <- coefs[1, 2]
      se_slope     <- coefs[2, 2]
      term1 <- (se_intercept * 1.96 / intercept_val) ^ 2
      term2 <- (se_slope * 1.96 / slope_val) ^ 2
      CI_val <- sqrt(term1 + term2) * abs(intercept_val / slope_val)
    }

    return(list(
      pvalue = p_val,
      pseudoR2 = pseudoR2_val,
      estimate = est_val,
      switch_at_time = switch_time_val,
      prd_quality = quality_val,
      CI = CI_val
    ))
  }

  # Run parallelized loop
  gene_indices <- 1:ncol(binary_t)

  if (.Platform$OS.type == "windows") {
    results_list <- pbapply::pblapply(gene_indices, fit_single_gene, cl = ncores)
  } else {
    results_list <- pbmcapply::pbmclapply(gene_indices, fit_single_gene, mc.cores = ncores)
  }

  message("Finished fitting models. Processing results...")

  # SAFETY CHECK: Ensure results_list is not NULL
  if (is.null(results_list)) {
    stop("Parallel processing failed. results_list is NULL.")
  }

  # Define the fail structure
  fail_struct <- list(
    pvalue = NA_real_,
    pseudoR2 = NA_real_,
    estimate = NA_real_,
    switch_at_time = NA_real_,
    prd_quality = 0,
    CI = NA_real_
  )


  # Unpack results
  # We need to extract these into separate vectors for further processing.
  # data.table::rbindlist could be faster?

  # model fitting can produce wanrings
  # if there are any warnings they can be appended as a second element in results_list
  # if results_list$warnings exists, remove it before unpacking results
  if (is.list(results_list) && "warnings" %in% names(results_list)) {
    results_list <- results_list$value
  }

  # Use lapply + unlist instead of sapply to ensure we get a vector, not a matrix/list if simplification behaves oddly
  pvalues        <- sapply(results_list, function(x) x$pvalue)
  pseudoR2s      <- sapply(results_list, function(x) x$pseudoR2)
  estimates      <- sapply(results_list, function(x) x$estimate)
  switch_at_time <- sapply(results_list, function(x) x$switch_at_time)
  prd_quality    <- sapply(results_list, function(x) x$prd_quality)
  CI             <- sapply(results_list, function(x) x$CI)

  # --- POST-PROCESSING ---
  # Create Dataframe of Results
  # We use 'colnames(binary_t)' because our data was transposed
  result_switch <- data.frame(
    switch_at_time = switch_at_time,
    CI             = CI,
    pvalues        = pvalues,
    pseudoR2s      = pseudoR2s,
    estimates      = estimates,
    prd_quality    = prd_quality,
    row.names      = colnames(binary_t) # Ensures gene names match results
  )

  # Add Derived Columns (Vectorized Operations)
  # Use ifelse for "up/down" logic—it's cleaner and faster
  result_switch$direction <- ifelse(result_switch$estimates > 0, "up", "down")
  result_switch$FDR       <- p.adjust(result_switch$pvalues, method = "BH")

  # Use our pre-calculated t_max and t_min constants
  steptime <- (t_max - t_min) / 100
  result_switch$switch_at_timeidx <- round((result_switch$switch_at_time - t_min) / steptime)

  # Quality Filter
  # Identify genes with high FDR.
  # Use which() because it handles NAs gracefully.
  # (If FDR is NA, which() simply ignores it, preventing the 'if' crash)
  high_fdr_indices <- which(result_switch$FDR > sig_FDR)
  if (length(high_fdr_indices) > 0) {
    result_switch$prd_quality[high_fdr_indices] <- 0
  }


  # --- INTEGRATION ---
  # Instead of a slow 'merge', we assign values directly to the SCE object.
  # This preserves row order and handles the fact that we filtered genes earlier.

  # Extract rowData to a temp variable (modifying S4 objects repeatedly is slow)
  r_data <- rowData(sce)

  # Initialize columns with NA (or 0) for ALL genes
  r_data$switch_at_time     <- NA_real_
  r_data$CI                 <- NA_real_
  r_data$pvalues            <- NA_real_
  r_data$pseudoR2s          <- NA_real_
  r_data$estimates          <- NA_real_
  r_data$prd_quality        <- 0  # Default quality is 0 (low) until proven otherwise
  r_data$direction          <- NA_character_
  r_data$FDR                <- NA_real_
  r_data$switch_at_timeidx  <- NA_real_


  # Inject results into the processed genes only
  r_data$switch_at_time    <- result_switch$switch_at_time
  r_data$CI                <- result_switch$CI
  r_data$pvalues           <- result_switch$pvalues
  r_data$pseudoR2s         <- result_switch$pseudoR2s
  r_data$estimates         <- result_switch$estimates
  r_data$prd_quality       <- result_switch$prd_quality
  r_data$direction         <- result_switch$direction
  r_data$FDR               <- result_switch$FDR
  r_data$switch_at_timeidx <- result_switch$switch_at_timeidx

  # Save back to SCE
  rowData(sce) <- r_data

  return(sce)
}
