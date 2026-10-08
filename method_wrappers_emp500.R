#!/usr/bin/env Rscript
# method_wrappers_emp500.R
# Wrapper functions for 15 differential abundance methods on EMP500 data.
# Formula: ~ empo_3 + env_biome
# empo_3 has 14 levels, env_biome has 35 levels.
#
# Each wrapper takes: counts (taxa x samples), metadata (samples x vars)
# Returns: data.frame with columns: taxon, p_value, q_value

# ============================================================
# Helper functions
# ============================================================

set_single_threaded <- function() {
  # 1. C/C++ System Environment Variables
  # Set these before libraries load to catch OpenMP and BLAS backends
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  
  # 2. BiocParallel
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam())
  }
  
  # 3. RcppParallel
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    RcppParallel::setThreadOptions(numThreads = 1)
  } 
  
  # 4. data.table
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(1)
  }
  
  # 5. BLAS / LAPACK / OpenMP
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  } else {
    warning("Package 'RhpcBLASctl' is not installed. BLAS multithreading might still occur.")
  }
}

# EMP500 formula (with ~ for model.matrix-based methods)
get_formula_str <- function() {
  "~ empo_3 + env_biome"
}

# Formula without ~ (for ANCOMBC, ANCOMBC2)
get_formula_str_no_tilde <- function() {
  "empo_3 + env_biome"
}

get_adj_vars <- function() {
  c("env_biome")
}

get_primary_var <- function() {
  "empo_3"
}

# Preprocess EMP500 data: remove zero-count samples, group rare factor levels,
# resolve collinearity between empo_3 and env_biome
preprocess_emp500 <- function(counts, metadata, min_level_size = 5) {
  # Remove samples with zero total counts
  lib_sizes <- colSums(counts)
  if (any(lib_sizes == 0)) {
    keep <- lib_sizes > 0
    counts <- counts[, keep, drop = FALSE]
    metadata <- metadata[keep, , drop = FALSE]
  }

  # Group rare empo_3 levels as "Other"
  primary_var <- get_primary_var()
  tbl <- table(metadata[[primary_var]])
  rare_levels <- names(tbl)[tbl < min_level_size]
  if (length(rare_levels) > 0) {
    vals <- as.character(metadata[[primary_var]])
    vals[vals %in% rare_levels] <- "Other"
    metadata[[primary_var]] <- factor(vals)
  }

  # Group rare env_biome levels as "Other"
  tbl2 <- table(metadata$env_biome)
  rare_biomes <- names(tbl2)[tbl2 < min_level_size]
  if (length(rare_biomes) > 0) {
    vals <- as.character(metadata$env_biome)
    vals[vals %in% rare_biomes] <- "Other"
    metadata$env_biome <- factor(vals)
  }

  # Resolve collinearity: iteratively merge env_biome levels that cause
  # non-full-rank design matrices into "Other"
  for (iter in 1:10) {
    design <- model.matrix(~ empo_3 + env_biome, data = as.data.frame(metadata))
    qr_obj <- qr(design)
    if (qr_obj$rank == ncol(design)) break  # full rank

    # Find non-estimable columns (env_biome only, keep empo_3)
    pivot <- qr_obj$pivot
    non_est <- pivot[(qr_obj$rank + 1):ncol(design)]
    non_est_cols <- colnames(design)[non_est]
    # Filter to env_biome columns only
    biome_cols <- non_est_cols[grepl("^env_biome", non_est_cols)]
    if (length(biome_cols) == 0) break

    # Extract level names from column names (e.g., "env_biomeXYZ" -> "XYZ")
    biome_levels_to_merge <- sub("^env_biome", "", biome_cols)
    vals <- as.character(metadata$env_biome)
    vals[vals %in% biome_levels_to_merge] <- "Other"
    metadata$env_biome <- factor(vals)
  }

  list(counts = counts, metadata = metadata)
}

# ============================================================
# 1. ORM (Ordinal Regression Model) - rms::orm
# ============================================================
run_orm_method <- function(counts, metadata) {
  set_single_threaded()
  library(rms)
  library(purrr)
  library(dplyr)

  rel_abn <- sweep(counts, 2, colSums(counts), "/")
  rel_abn[is.na(rel_abn)] <- 0

  run_orm <- function(abundance, metadata, formula) {
    tryCatch({
    mm <- model.matrix(formula, metadata) |>
      cbind.data.frame(abundance)
    mm[["(Intercept)"]] <- NULL
    inds <- 1:(ncol(mm) - 1)
    vars <- colnames(mm)[inds]
    fit_1 <- tryCatch(orm(abundance ~ ., data = mm),
                      error = function(e) return(NULL))
    if (is.null(fit_1)) return(NULL)
    score_1 <- fit_1$stats["Score"]
    res <- data.frame(estimate = fit_1$coefficients[vars], p_value = NA)
    if (length(inds) == 1) {
      res$p_value <- fit_1$stats["Score P"]
    } else {
      for (i in inds) {
        fit_0 <- tryCatch(orm(abundance ~ ., data = mm[, -i]),
                          error = function(e) return(NULL))
        if (!is.null(fit_0)) {
          score_0 <- fit_0$stats["Score"]
          res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
        }
      }
    }
    res[["variable"]] <- rownames(res)
    rownames(res) <- NULL
    res
    }, error = function(e) NULL)
  }

  formula <- as.formula(get_formula_str())
  meta_df <- as.data.frame(metadata)
  primary_var <- get_primary_var()

  res <- as.data.frame(t(rel_abn)) |>
    map(~ run_orm(., metadata = meta_df, formula = formula)) |>
    compact() |>
    bind_rows(.id = "taxon") |>
    filter(grepl(paste0("^", primary_var), variable)) |>
    group_by(taxon) |>
    summarise(p_value = min(p.adjust(p_value, method = "bonferroni")), .groups = "drop") |>
    mutate(q_value = p.adjust(p_value, method = "BH")) |>
    select(taxon, p_value, q_value)

  return(res)
}

# ============================================================
# 2. MaAsLin2
# ============================================================
run_maaslin2 <- function(counts, metadata) {
  set_single_threaded()
  library(Maaslin2)

  tmp_dir <- tempfile("maaslin2_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  counts_df <- as.data.frame(counts)
  meta_df <- as.data.frame(metadata)
  fixed_effects <- colnames(meta_df)
  primary_var <- get_primary_var()

  tryCatch({
    Maaslin2(
      input_data = counts_df,
      input_metadata = meta_df,
      output = tmp_dir,
      fixed_effects = fixed_effects,
      random_effects = NULL,
      reference = paste0(primary_var, ",", levels(meta_df[[primary_var]])[1]),
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      correction = "BH",
      min_abundance = 0,
      min_prevalence = 0,
      cores = 1
    )
  }, error = function(e) {
    if (!file.exists(file.path(tmp_dir, "all_results.tsv")))
      stop(paste("MaAsLin2 failed:", conditionMessage(e)))
  })

  res_df <- read.csv(file.path(tmp_dir, "all_results.tsv"), sep = "\t")
  res_df <- res_df[res_df$metadata == primary_var, ]
  res_df <- res_df %>%
    group_by(feature) %>%
    summarise(p_value = min(p.adjust(pval, method = "bonferroni")), .groups = "drop") %>%
    mutate(q_value = p.adjust(p_value, method = "BH"))
  res_df <- data.frame(taxon = res_df$feature, p_value = res_df$p_value, q_value = res_df$q_value)

  return(res_df)
}

# ============================================================
# 3. MaAsLin3
# ============================================================
run_maaslin3 <- function(counts, metadata) {
  set_single_threaded()
  library(maaslin3)
  library(dplyr)

  tmp_dir <- tempfile("maaslin3_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  counts_df <- as.data.frame(t(as.matrix(counts)))
  meta_df <- as.data.frame(metadata)
  formula <- get_formula_str()
  primary_var <- get_primary_var()

  res <- maaslin3(
    input_data = counts_df,
    input_metadata = meta_df,
    output = tmp_dir,
    formula = formula,
    reference = paste0(primary_var, ",", levels(meta_df[[primary_var]])[1]),
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    min_abundance = 0,
    min_prevalence = 0,
    standardize = FALSE,
    augment = FALSE,
    max_significance = 0.05,
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    save_models = FALSE,
    cores = 1,
    verbosity = "WARN"
  )

  res_file <- file.path(tmp_dir, "all_results.tsv")
  if (file.exists(res_file)) {
    res_df <- read.csv(res_file, sep = "\t")
    res_df <- res_df[res_df$metadata == primary_var & res_df$model == "abundance", ]
    res_df <- res_df %>%
      group_by(feature) %>%
      summarise(p_value = min(p.adjust(pval_individual, method = "bonferroni")), .groups = "drop") %>%
      mutate(q_value = p.adjust(p_value, method = "BH"))
    res_df <- data.frame(taxon = res_df$feature, p_value = res_df$p_value, q_value = res_df$q_value)
  } else {
    res_df <- data.frame(taxon = rownames(counts), p_value = NA, q_value = NA)
  }

  return(res_df)
}

# ============================================================
# 4. LinDA
# ============================================================
run_linda <- function(counts, metadata) {
  set_single_threaded()
  library(LinDA)

  meta_df <- as.data.frame(metadata)
  formula <- get_formula_str()
  primary_var <- get_primary_var()

  res <- linda(
    as.matrix(counts),
    meta_df,
    formula = formula,
    alpha = 0.05,
    prev.cut = 0,
    lib.cut = 0,
    adaptive = FALSE,
    imputation = FALSE,
    winsor.quan = NULL
  )

  coef_names <- names(res$output)
  primary_coefs <- grep(paste0("^", primary_var), coef_names, value = TRUE)

  if (length(primary_coefs) > 0) {
    all_pvals <- do.call(rbind, lapply(primary_coefs, function(cn) {
      tab <- res$output[[cn]]
      data.frame(taxon = rownames(tab), pval = tab$pval, stringsAsFactors = FALSE)
    }))
    res_df <- aggregate(pval ~ taxon, data = all_pvals,
                        FUN = function(x) min(p.adjust(x, method = "bonferroni")))
    colnames(res_df) <- c("taxon", "p_value")
    res_df$q_value <- p.adjust(res_df$p_value, method = "BH")
  } else {
    res_df <- data.frame(taxon = rownames(counts), p_value = NA, q_value = NA)
  }

  return(res_df)
}

# ============================================================
# 5. corncob
# ============================================================
run_corncob <- function(counts, metadata) {
  set_single_threaded()
  library(corncob)
  library(phyloseq)

  otu_tab <- otu_table(counts, taxa_are_rows = TRUE)
  tax_tab <- tax_table(matrix("unknown", nrow = nrow(counts), ncol = 1,
                               dimnames = list(rownames(counts), "Species")))
  sample_data <- sample_data(metadata)

  ps <- phyloseq(otu_tab, tax_tab, sample_data)

  formula <- as.formula(get_formula_str())
  phi_formula <- as.formula("~ 1")
  formula_null <- as.formula("~ env_biome")

  res <- tryCatch(
    differentialTest(
      formula = formula,
      phi.formula = phi_formula,
      phi.formula_null = phi_formula,
      formula_null = formula_null,
      data = ps,
      test = "Wald",
      boot = FALSE,
      fdr = "BH"
    ),
    error = function(e) {
      # Retry with simpler formula (empo_3 only, no env_biome adj)
      tryCatch(
        differentialTest(
          formula = formula,
          phi.formula = phi_formula,
          phi.formula_null = phi_formula,
          formula_null = ~ 1,
          data = ps,
          test = "Wald",
          boot = FALSE,
          fdr = "BH"
        ),
        error = function(e2) {
          # Last resort: empo_3 only as both formula and null
          differentialTest(
            formula = ~ empo_3,
            phi.formula = phi_formula,
            phi.formula_null = phi_formula,
            formula_null = ~ 1,
            data = ps,
            test = "Wald",
            boot = FALSE,
            fdr = "BH"
          )
        }
      )
    }
  )

  res_df <- data.frame(
    taxon = names(res$p),
    p_value = res$p,
    q_value = p.adjust(res$p, method = "BH")
  )

  return(res_df)
}

# ============================================================
# 6. limma-voom
# ============================================================
run_limma_voom <- function(counts, metadata) {
  set_single_threaded()
  library(limma)
  library(edgeR)

  meta_df <- as.data.frame(metadata)
  design <- model.matrix(as.formula(get_formula_str()), data = meta_df)

  dge <- DGEList(counts = as.matrix(counts))
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  primary_var <- get_primary_var()
  primary_cols <- grep(paste0("^", primary_var), colnames(design))

  tt <- topTable(fit, coef = primary_cols, number = Inf, sort.by = "none")

  res_df <- data.frame(
    taxon = rownames(tt),
    p_value = tt$P.Value,
    q_value = tt$adj.P.Val
  )

  return(res_df)
}

# ============================================================
# 7. ALDEx2
# ============================================================
run_aldex2 <- function(counts, metadata) {
  set_single_threaded()
  library(ALDEx2)

  counts_num <- matrix(as.double(as.matrix(counts)),
                       nrow = nrow(counts), ncol = ncol(counts))
  rownames(counts_num) <- rownames(counts)
  colnames(counts_num) <- colnames(counts)

  primary_var <- get_primary_var()
  primary_vals <- metadata[[primary_var]]
  most_common <- names(sort(table(primary_vals), decreasing = TRUE))[1]
  conditions <- ifelse(as.character(primary_vals) == most_common, "A", "B")
  conditions <- as.character(conditions)

  res <- aldex(
    reads = counts_num,
    conditions = conditions,
    mc.samples = 128,
    test = "t",
    effect = FALSE,
    verbose = FALSE
  )

  res_df <- data.frame(
    taxon = rownames(res),
    p_value = res$we.eBH,
    q_value = res$we.eBH
  )

  return(res_df)
}

# ============================================================
# 8. ZicoSeq
# ============================================================
run_zicoseq <- function(counts, metadata) {
  set_single_threaded()
  library(GUniFrac)

  meta_df <- as.data.frame(metadata)
  primary_var <- get_primary_var()

  most_common <- names(sort(table(meta_df[[primary_var]]), decreasing = TRUE))[1]
  meta_df$empo_3_binary <- ifelse(as.character(meta_df[[primary_var]]) == most_common, "A", "B")

  res <- ZicoSeq(
    meta.dat = meta_df,
    feature.dat = as.matrix(counts),
    grp.name = "empo_3_binary",
    adj.name = "env_biome",
    feature.dat.type = "count",
    prev.filter = 0,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    min.prop = 0,
    is.winsor = FALSE,
    outlier.pct = 0.03,
    winsor.end = "top",
    is.post.sample = TRUE,
    post.sample.no = 25,
    link.func = list(function(x) sign(x) * (abs(x))^0.5),
    stats.combine.func = max,
    perm.no = 99,
    strata = NULL,
    ref.pct = 0.5,
    stage.no = 6,
    excl.pct = 0.2,
    is.fwer = FALSE,
    verbose = FALSE
  )

  pvals <- res$p.raw
  res_df <- data.frame(
    taxon = names(pvals),
    p_value = pvals,
    q_value = p.adjust(pvals, method = "BH")
  )

  return(res_df)
}

# ============================================================
# 9. LDM
# fixed sorce code in LDM_fun.R: rownames(beta) = make.unique(beta.name)
# ============================================================
run_ldm <- function(counts, metadata) {
  set_single_threaded()
  library(LDM)
#  meta_df <- as.data.frame(metadata)
#  otu_tab <- t(as.matrix(counts))
#
#  # LDM has issues with factor level names containing spaces/parentheses
#  # Save original levels and temporarily replace with simple names
#  orig_empo_levels <- levels(meta_df$empo_3)
#  levels(meta_df$empo_3) <- paste0("E", seq_along(orig_empo_levels))
#
#  .otu_tab_existed <- exists(".ldm_otu_tab", envir = globalenv())
#  if (.otu_tab_existed) .old_otu_tab <- get(".ldm_otu_tab", envir = globalenv())
#  assign(".ldm_otu_tab", otu_tab, envir = globalenv())
#
#  formula_str <- ".ldm_otu_tab ~ empo_3"
# Temporarily silence the matrixStats ties.method warning
    options(matrixStats.options = list(ties.method.warn = FALSE))

    meta_df <- as.data.frame(metadata)
    otu_tab <- t(as.matrix(counts))

    # 1. Clean the factor: make sure it's a factor and drop any empty/unused levels
    meta_df$empo_3 <- as.factor(meta_df$empo_3)
    meta_df$empo_3 <- droplevels(meta_df$empo_3) 

    # LDM has issues with factor level names containing spaces/parentheses
    # Save original levels and temporarily replace with simple names
    orig_empo_levels <- levels(meta_df$empo_3)
    levels(meta_df$empo_3) <- paste0("Grp_", seq_along(orig_empo_levels))

    # 2. Restore your exact environment logic so the downstream script doesn't break
    .otu_tab_existed <- exists(".ldm_otu_tab", envir = globalenv())
    if (.otu_tab_existed) {
        .old_otu_tab <- get(".ldm_otu_tab", envir = globalenv())
    }

    # Assign the new one
    assign(".ldm_otu_tab", otu_tab, envir = globalenv())

    formula_str <- ".ldm_otu_tab ~ empo_3"

  tryCatch({
    res <- ldm(
      formula = as.formula(formula_str),
      data = meta_df,
      n.perm.max = 1000,
      n.cores = 1,
      verbose = FALSE
    )
  }, error = function(e) {
    stop(paste("LDM failed:", conditionMessage(e)))
  }, finally = {
    if (.otu_tab_existed) assign(".ldm_otu_tab", .old_otu_tab, envir = globalenv())
    else rm(".ldm_otu_tab", envir = globalenv())
  })

  if (!is.null(res$p.otu.tran)) {
    pvals <- res$p.otu.tran[1, ]
    qvals <- if (!is.null(res$q.otu.tran)) res$q.otu.tran[1, ] else p.adjust(pvals, method = "BH")
  } else if (!is.null(res$p.otu.freq)) {
    pvals <- res$p.otu.freq[1, ]
    qvals <- if (!is.null(res$q.otu.freq)) res$q.otu.freq[1, ] else p.adjust(pvals, method = "BH")
  } else {
    pvals <- rep(NA, ncol(otu_tab))
    qvals <- rep(NA, ncol(otu_tab))
  }

  res_df <- data.frame(
    taxon = colnames(otu_tab),
    p_value = pvals,
    q_value = qvals
  )

  return(res_df)
}

# ============================================================
# 10. metagenomeSeq
# ============================================================
run_metagenomeseq <- function(counts, metadata) {
  set_single_threaded()
  library(metagenomeSeq)
  library(Biobase)

  meta_df <- as.data.frame(metadata)
  pd <- AnnotatedDataFrame(meta_df)
  obj <- newMRexperiment(as.matrix(counts), phenoData = pd)

  # cumNorm can fail on sparse subsets; use manual normalization as fallback
  obj <- tryCatch(
    cumNorm(obj),
    error = function(e) {
      # Manual CSS normalization: use median library size as scaling factor
      nf <- colSums(as.matrix(counts))
      nf <- nf / median(nf[nf > 0])
      nf[!is.finite(nf)] <- 1
      nf[nf <= 0] <- 1
      normFactors(obj) <- log2(nf + 1)
      obj
    }
  )

  design <- model.matrix(as.formula(get_formula_str()), data = meta_df)
  primary_var <- get_primary_var()

  # fitFeatureModel often fails with complex designs; use fitZig directly
  fit <- fitZig(obj, mod = design, useCSSoffset = TRUE,
                control = zigControl(verbose = FALSE))

  # Extract per-coefficient p-values using MRfulltable
  mr_full <- tryCatch(MRfulltable(fit, number = Inf), error = function(e) NULL)

  if (!is.null(mr_full)) {
    pval_cols <- grep("\\+pvalue$", colnames(mr_full), value = TRUE)
    primary_pcols <- pval_cols[grep(primary_var, pval_cols)]
    if (length(primary_pcols) > 0) {
      all_pvals <- mr_full[, primary_pcols, drop = FALSE]
      combined_p <- apply(all_pvals, 1, function(x) min(p.adjust(x, method = "bonferroni")))
      res_df <- data.frame(
        taxon = rownames(mr_full),
        p_value = combined_p,
        q_value = p.adjust(combined_p, method = "BH")
      )
      res_df <- res_df[!is.na(res_df$p_value), ]
      return(res_df)
    }
  }

  # Fallback: use MRtable (overall p-value)
  mr <- MRtable(fit, number = Inf)
  res_df <- data.frame(
    taxon = rownames(mr),
    p_value = mr$pvalues,
    q_value = mr$adjPvalues
  )
  res_df <- res_df[!is.na(res_df$p_value), ]
  return(res_df)
}

# ============================================================
# 11. edgeR
# ============================================================
run_edger <- function(counts, metadata) {
  set_single_threaded()
  library(edgeR)

  meta_df <- as.data.frame(metadata)
  design <- model.matrix(as.formula(get_formula_str()), data = meta_df)

  dge <- DGEList(counts = as.matrix(counts))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  # estimateDisp fails on extremely sparse microbiome data (97% zeros)
  # due to locfit crashing in the trend dispersion estimation.
  # Use a fixed common dispersion (0.1) instead — a reasonable default
  # for microbial count data.
  dge$common.dispersion <- 0.1

  fit <- glmFit(dge, design, dispersion = 0.1)
  primary_var <- get_primary_var()
  primary_cols <- grep(paste0("^", primary_var), colnames(design))
  lrt <- glmLRT(fit, coef = primary_cols)

  tt <- topTags(lrt, n = Inf, sort.by = "none")$table

  res_df <- data.frame(
    taxon = rownames(tt),
    p_value = tt$PValue,
    q_value = tt$FDR
  )

  return(res_df)
}

# ============================================================
# 12. DESeq2
# ============================================================
run_deseq2 <- function(counts, metadata) {
  set_single_threaded()
  library(DESeq2)

  meta_df <- as.data.frame(metadata)

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts),
    colData = meta_df,
    design = as.formula(get_formula_str())
  )

  dds <- estimateSizeFactors(dds, type = "poscounts")

  # LRT: full model ~ empo_3 + env_biome vs reduced ~ env_biome
  dds <- DESeq(dds, test = "LRT",
               reduced = as.formula("~ env_biome"),
               parallel = FALSE)

  res <- results(dds, alpha = 0.05)

  res_df <- data.frame(
    taxon = rownames(res),
    p_value = res$pvalue,
    q_value = res$padj
  )
  res_df <- res_df[!is.na(res_df$p_value), ]

  return(res_df)
}

# ============================================================
# 13. ANCOMBC
# ============================================================
run_ancombc <- function(counts, metadata) {
  set_single_threaded()
  library(ANCOMBC)
  library(phyloseq)
  library(microbiome)

  # Create phyloseq object
  otu_tab <- otu_table(as.matrix(counts), taxa_are_rows = TRUE)
  sample_data <- sample_data(as.data.frame(metadata))
  ps <- phyloseq(otu_tab, sample_data)

  formula_str <- get_formula_str_no_tilde()
  primary_var <- get_primary_var()

  res <- ancombc(
    data = ps,
    formula = formula_str,
    p_adj_method = "holm",
    prv_cut = 0.0,
    lib_cut = 0,
    struc_zero = FALSE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = FALSE,
    alpha = 0.05,
    n_cl = 1,
    global = FALSE,
    verbose = FALSE
  )

  # res$res is a list with elements: lfc, se, W, p_val, q_val, diff_abn
  # p_val is a data.frame: taxon, (Intercept), empo_3Level1, empo_3Level2, ..., env_biomeLevel1, ...
  p_val_df <- res$res$p_val
  q_val_df <- res$res$q_val

  # Find columns matching the primary variable (empo_3)
  primary_cols <- grep(paste0("^", primary_var), colnames(p_val_df), value = TRUE)

  if (length(primary_cols) > 0) {
    # Combine p-values across empo_3 coefficients per taxon (Bonferroni)
    pvals_mat <- as.matrix(p_val_df[, primary_cols, drop = FALSE])
    combined_p <- apply(pvals_mat, 1, function(x) min(p.adjust(x, method = "bonferroni")))
    res_df <- data.frame(
      taxon = p_val_df$taxon,
      p_value = combined_p,
      q_value = p.adjust(combined_p, method = "BH")
    )
  } else {
    res_df <- data.frame(taxon = rownames(counts), p_value = NA, q_value = NA)
  }

  return(res_df)
}

# ============================================================
# 14. ANCOMBC2
# ============================================================
run_ancombc2 <- function(counts, metadata) {
  set_single_threaded()
  library(ANCOMBC)
  library(phyloseq)
  library(microbiome)

  # Create phyloseq object
  otu_tab <- otu_table(as.matrix(counts), taxa_are_rows = TRUE)
  sample_data <- sample_data(as.data.frame(metadata))
  ps <- phyloseq(otu_tab, sample_data)

  formula_str <- get_formula_str_no_tilde()
  primary_var <- get_primary_var()

  res <- ancombc2(
    data = ps,
    fix_formula = formula_str,
    rand_formula = NULL,
    p_adj_method = "holm",
    pseudo = 0,
    pseudo_sens = TRUE,
    prv_cut = 0.0,
    lib_cut = 0,
    s0_perc = 0.02,
    struc_zero = FALSE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    verbose = FALSE
  )

  # res$res is a data.frame with columns: taxon, lfc_*, se_*, W_*, p_*, q_*, diff_*, passed_ss_*
  res_df_main <- res$res

  # Find p-value columns matching the primary variable
  pval_cols <- grep(paste0("^p_", primary_var), colnames(res_df_main), value = TRUE)

  if (length(pval_cols) > 0) {
    pvals_mat <- as.matrix(res_df_main[, pval_cols, drop = FALSE])
    combined_p <- apply(pvals_mat, 1, function(x) min(p.adjust(x, method = "bonferroni")))
    res_df <- data.frame(
      taxon = res_df_main$taxon,
      p_value = combined_p,
      q_value = p.adjust(combined_p, method = "BH")
    )
  } else {
    res_df <- data.frame(taxon = rownames(counts), p_value = NA, q_value = NA)
  }

  return(res_df)
}

# ============================================================
# Method registry
# ============================================================
get_method_functions <- function() {
  list(
    orm = run_orm_method,
    maaslin2 = run_maaslin2,
    maaslin3 = run_maaslin3,
    linda = run_linda,
    corncob = run_corncob,
    limma_voom = run_limma_voom,
    aldex2 = run_aldex2,
    zicoseq = run_zicoseq,
    ldm = run_ldm,
    metagenomeseq = run_metagenomeseq,
    edger = run_edger,
    deseq2 = run_deseq2,
    ancombc = run_ancombc,
    ancombc2 = run_ancombc2
  )
}

get_method_packages <- function() {
  list(
    orm = "rms",
    maaslin2 = "Maaslin2",
    maaslin3 = "maaslin3",
    linda = "LinDA",
    corncob = "corncob",
    limma_voom = "limma",
    aldex2 = "ALDEx2",
    zicoseq = "GUniFrac",
    ldm = "LDM",
    metagenomeseq = "metagenomeSeq",
    edger = "edgeR",
    deseq2 = "DESeq2",
    ancombc = "ANCOMBC",
    ancombc2 = "ANCOMBC"
  )
}
