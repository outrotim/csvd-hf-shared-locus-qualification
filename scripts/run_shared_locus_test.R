#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
usage <- "usage: run_shared_locus_test.R exposure.tsv outcome.tsv output_prefix abf|susie [ld_matrix.rds]"
if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
  cat(usage, "\n")
  quit(status = 0L)
}
if (length(args) < 4L) stop(usage)

exposure_path <- args[[1L]]
outcome_path <- args[[2L]]
output_prefix <- args[[3L]]
mode <- args[[4L]]
ld_path <- if (length(args) >= 5L) args[[5L]] else NA_character_

required <- c("snp", "beta", "se", "maf", "position", "n")
read_trait <- function(path, label) {
  x <- fread(path)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ","))
  x <- x[
    is.finite(beta) & is.finite(se) & se > 0 &
      is.finite(maf) & maf > 0 & maf <= 0.5 &
      is.finite(position) & is.finite(n) & n > 0
  ]
  if (anyDuplicated(x$snp)) stop(label, " contains duplicate SNP identifiers")
  x
}

exposure <- read_trait(exposure_path, "exposure")
outcome <- read_trait(outcome_path, "outcome")
shared <- merge(exposure, outcome, by = "snp", suffixes = c("_exposure", "_outcome"))
setorder(shared, position_exposure, snp)
if (nrow(shared) < 100L) stop("fewer than 100 finite shared variants")
if (max(abs(shared$position_exposure - shared$position_outcome)) > 0) {
  stop("exposure/outcome position mismatch after build harmonization")
}

make_exposure <- function(rows) {
  list(
    beta = rows$beta_exposure,
    varbeta = rows$se_exposure^2,
    snp = rows$snp,
    position = rows$position_exposure,
    MAF = rows$maf_exposure,
    type = "quant",
    N = median(rows$n_exposure),
    sdY = 1
  )
}

make_outcome <- function(rows, case_fraction = 0.07) {
  list(
    beta = rows$beta_outcome,
    varbeta = rows$se_outcome^2,
    snp = rows$snp,
    position = rows$position_outcome,
    MAF = rows$maf_outcome,
    type = "cc",
    N = median(rows$n_outcome),
    s = case_fraction
  )
}

dataset1 <- make_exposure(shared)
dataset2 <- make_outcome(shared, 0.07)
check_dataset(dataset1)
check_dataset(dataset2)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
summary_rows <- list()
pair_table <- data.table()

if (mode == "abf") {
  for (fraction in c(0.06, 0.07, 0.08)) {
    for (p12 in c(5e-6, 1e-5, 5e-5)) {
      fit <- coloc.abf(
        dataset1 = dataset1,
        dataset2 = make_outcome(shared, fraction),
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = p12
      )
      s <- fit$summary
      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        mode = "abf",
        case_fraction = fraction,
        p12 = p12,
        nsnps = as.numeric(s[["nsnps"]]),
        pp0 = as.numeric(s[["PP.H0.abf"]]),
        pp1 = as.numeric(s[["PP.H1.abf"]]),
        pp2 = as.numeric(s[["PP.H2.abf"]]),
        pp3 = as.numeric(s[["PP.H3.abf"]]),
        pp4 = as.numeric(s[["PP.H4.abf"]]),
        primary = fraction == 0.07 && p12 == 1e-5
      )
    }
  }
} else if (mode == "susie") {
  if (is.na(ld_path) || !file.exists(ld_path)) stop("SuSiE mode requires an LD matrix RDS")
  ld <- readRDS(ld_path)
  if (!is.matrix(ld) || is.null(rownames(ld)) || is.null(colnames(ld))) {
    stop("LD must be a named square matrix")
  }
  if (!identical(rownames(ld), colnames(ld))) stop("LD row/column names differ")
  if (!all(shared$snp %in% rownames(ld))) stop("shared variants absent from LD matrix")
  ld <- ld[shared$snp, shared$snp, drop = FALSE]
  if (max(abs(ld - t(ld))) > 1e-8 || max(abs(diag(ld) - 1)) > 1e-8) {
    stop("LD matrix is not a valid signed correlation matrix")
  }
  dataset1$LD <- ld
  fit1 <- runsusie(dataset1, coverage = 0.95, maxit = 1000)
  pair_tables <- list()
  for (fraction in c(0.06, 0.07, 0.08)) {
    outcome_dataset <- make_outcome(shared, fraction)
    outcome_dataset$LD <- ld
    fit2 <- runsusie(outcome_dataset, coverage = 0.95, maxit = 1000)
    for (p12 in c(5e-6, 1e-5, 5e-5)) {
      fit <- coloc.susie(fit1, fit2, p1 = 1e-4, p2 = 1e-4, p12 = p12)
      table <- as.data.table(fit$summary)
      if (!nrow(table)) next
      old <- intersect(c("PP.H3.abf", "PP.H4.abf"), names(table))
      new <- c("PP.H3.abf" = "pp3", "PP.H4.abf" = "pp4")[old]
      setnames(table, old = old, new = unname(new))
      table[, `:=`(
        case_fraction = fraction,
        p12 = p12,
        primary = fraction == 0.07 && p12 == 1e-5
      )]
      pair_tables[[length(pair_tables) + 1L]] <- table
    }
  }
  if (!length(pair_tables)) {
    pair_table <- data.table(
      mode = "susie",
      case_fraction = 0.07,
      p12 = 1e-5,
      nsnps = nrow(shared),
      pp3 = NA_real_,
      pp4 = NA_real_,
      primary = TRUE,
      status = "NO_CREDIBLE_SET_PAIR"
    )
  } else {
    pair_table <- rbindlist(pair_tables, fill = TRUE)
    pair_table[, status := "CREDIBLE_SET_PAIR_ESTIMATED"]
  }
} else {
  stop("mode must be 'abf' or 'susie'")
}

summary_table <- if (length(summary_rows)) rbindlist(summary_rows) else pair_table
if (!all(c("pp3", "pp4") %in% names(summary_table))) {
  stop("colocalization result lacks PP3/PP4")
}
summary_table[, gate_pass := primary & is.finite(pp4) & is.finite(pp3) & pp4 >= 0.80 & pp3 < 0.10]
fwrite(summary_table, paste0(output_prefix, ".summary.tsv"), sep = "\t")
fwrite(shared, paste0(output_prefix, ".shared_variants.tsv.gz"), sep = "\t")

no_pair <- mode == "susie" && all(summary_table$status == "NO_CREDIBLE_SET_PAIR")
audit <- list(
  state = if (no_pair) "FROZEN_COLOCALIZATION_COMPLETE_NO_CREDIBLE_SET_PAIR" else "FROZEN_COLOCALIZATION_COMPLETE",
  mode = mode,
  shared_variants = nrow(shared),
  primary_gate_pass = any(summary_table$gate_pass),
  credible_set_pair_rows = if (mode == "susie") sum(is.finite(summary_table$pp4)) else NA_integer_,
  negative_reason = if (no_pair) "NO_CREDIBLE_SET_PAIR" else NA_character_,
  p1 = 1e-4,
  p2 = 1e-4,
  primary_p12 = 1e-5,
  primary_case_fraction = 0.07,
  candidate_specific_selection_performed = FALSE
)
writeLines(toJSON(audit, auto_unbox = TRUE, pretty = TRUE), paste0(output_prefix, ".audit.json"))
print(audit)
