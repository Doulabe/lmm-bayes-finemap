# ==============================================================================
# run_all.R
# Pure-R orchestrator for the full benchmark suite. No bash, no GNU parallel.
# Uses parallel::mclapply (fork-based, macOS / Linux) to distribute work
# across local cores. Checkpointing per (cell, rep) RDS allows seamless resume.
#
# Usage from R:
#   source("sim/bench_full/run_all.R")
#   run_all(n_cores = 6, B = 10)                  # default: no BSLMM
#   run_all(n_cores = 8, B = 10, with_bslmm = TRUE)
#   run_all(n_cores = 6, B = 5, axes = c("01_n","05_mice"))   # subset
#   run_all(n_cores = 1)                          # serial
#
# Or from command line:
#   Rscript sim/bench_full/run_all.R              # defaults (n_cores=4, B=10)
#   Rscript sim/bench_full/run_all.R --n-cores 8 --B 5
#   Rscript sim/bench_full/run_all.R --axes 01_n,03_arch
#   Rscript sim/bench_full/run_all.R --with-bslmm --n-cores 8
#
# Design notes:
#   - Each (axis, cell, rep) is dispatched as one task to mclapply.
#   - already_done() checks the RDS path BEFORE doing any work, so resume
#     after interruption is automatic.
#   - On Windows, mclapply degrades to serial; pass n_cores=1 explicitly
#     and lean on the per-cell speed instead.
# ==============================================================================

suppressPackageStartupMessages({
  library(parallel)
})

# ------------------------------------------------------------------------------
# CRITICAL: limit BLAS threads to 1 per worker before forking.
# macOS Accelerate / OpenBLAS use multiple threads by default; with
# mclapply workers, each worker tries to use all cores -> heavy contention
# (we observed 30-60x slowdowns). Setting these env vars before R does any
# matrix algebra forces single-thread BLAS, so each mclapply worker gets a
# clean core. Effective only if the parent R process has not yet used BLAS.
# Override by setting LMM_BENCH_BLAS_THREADS in the environment if you know
# what you are doing.
# ------------------------------------------------------------------------------
.set_blas_threads <- function(nt = 1L) {
  Sys.setenv(OMP_NUM_THREADS         = nt,
             OPENBLAS_NUM_THREADS    = nt,
             MKL_NUM_THREADS         = nt,
             VECLIB_MAXIMUM_THREADS  = nt,
             BLIS_NUM_THREADS        = nt)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(nt)
    RhpcBLASctl::omp_set_num_threads(nt)
  }
}
.set_blas_threads(as.integer(Sys.getenv("LMM_BENCH_BLAS_THREADS", "1")))


# ------------------------------------------------------------------------------
# Working-directory-robust source(). Tries the candidate paths in order
# and uses whichever exists. Lets the user `source()` from anywhere:
#   - project root          (.../LMM_Gaussian_Bayesian/)
#   - sim/                  (.../LMM_Gaussian_Bayesian/sim/)
#   - sim/bench_full/       (.../LMM_Gaussian_Bayesian/sim/bench_full/)
# ------------------------------------------------------------------------------
.try_source_first <- function(candidates) {
  for (p in candidates) if (file.exists(p)) { source(p); return(p) }
  stop("Could not locate any of: ", paste(candidates, collapse = ", "),
       "\n  CWD = ", getwd(),
       "\n  Run from project root, sim/, or sim/bench_full/")
}
BENCH_DIR <- .try_source_first(c(
  "sim/bench_full/00_config.R",      # CWD = project root
  "bench_full/00_config.R",          # CWD = sim/
  "00_config.R"                       # CWD = sim/bench_full/
))
BENCH_DIR    <- normalizePath(dirname(BENCH_DIR))     # absolute path
PROJECT_ROOT <- normalizePath(file.path(BENCH_DIR, "..", ".."))


#=============================================================================
# Per-axis cell runners
# Each one takes the cell parameters + rep_id, generates data, runs methods,
# computes metrics, saves RDS. Returns invisible(NULL) on success.
# ============================================================================

# 01: vary n
run_cell_01n <- function(n, sg, rep_id, methods) {
  out_dir <- "results/bench_full/01_scaling_n"
  cell_tag <- sprintf("n%04d_sg%.1f", n, sg)
  if (already_done(out_dir, cell_tag, rep_id)) return("skip")
  seed <- 20260425L + 1000L * (rep_id - 1L) + as.integer(n) +
            round(1000 * sg)
  d <- gen_dataset(n = n, m = ANCHOR$m, rho = ANCHOR$rho,
                    K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- n; metrics$m <- ANCHOR$m; metrics$rho <- ANCHOR$rho
  metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  save_cell_rep(out_dir, cell_tag, rep_id,
                  list(metrics = metrics, truth = d$truth, results = res))
  "done"
}

# 02: vary rho
run_cell_02rho <- function(rho, sg, rep_id, methods) {
  out_dir <- "results/bench_full/02_scaling_rho"
  cell_tag <- sprintf("rho%.2f_sg%.1f", rho, sg)
  if (already_done(out_dir, cell_tag, rep_id)) return("skip")
  seed <- 20260425L + 1000L * (rep_id - 1L) +
            as.integer(round(1000 * rho)) + as.integer(round(100 * sg))
  d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = rho,
                    K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- ANCHOR$n; metrics$m <- ANCHOR$m; metrics$rho <- rho
  metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  save_cell_rep(out_dir, cell_tag, rep_id,
                  list(metrics = metrics, truth = d$truth, results = res))
  "done"
}

# 03: vary signal architecture
run_cell_03arch <- function(signal, sg, rep_id, methods) {
  out_dir <- "results/bench_full/03_arch"
  cell_tag <- sprintf("sig%s_sg%.1f", signal, sg)
  if (already_done(out_dir, cell_tag, rep_id)) return("skip")
  beta_true <- switch(signal,
                       weak    = c(0.4, 0.2, 0.2, 0.1, 0.1),
                       medium  = c(0.8, 0.4, 0.4, 0.2, 0.2),
                       strong  = c(1.6, 0.8, 0.8, 0.4, 0.4),
                       stop("Unknown signal: ", signal))
  seed <- 20260425L + 1000L * (rep_id - 1L) +
            switch(signal, weak = 1L, medium = 2L, strong = 3L) +
            as.integer(round(100 * sg))
  d <- gen_dataset(n = ANCHOR$n, m = ANCHOR$m, rho = ANCHOR$rho,
                    K_true = length(beta_true), beta_true = beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- ANCHOR$n; metrics$m <- ANCHOR$m; metrics$rho <- ANCHOR$rho
  metrics$signal <- signal; metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  save_cell_rep(out_dir, cell_tag, rep_id,
                  list(metrics = metrics, truth = d$truth,
                         results = res, beta_true = beta_true))
  "done"
}

# 04: vary m
run_cell_04m <- function(m, sg, rep_id, methods) {
  out_dir <- "results/bench_full/04_scaling_m"
  cell_tag <- sprintf("m%05d_sg%.1f", m, sg)
  if (already_done(out_dir, cell_tag, rep_id)) return("skip")
  seed <- 20260425L + 1000L * (rep_id - 1L) +
            as.integer(m) + as.integer(round(100 * sg))
  d <- gen_dataset(n = ANCHOR$n, m = m, rho = ANCHOR$rho,
                    K_true = ANCHOR$K_true, beta_true = ANCHOR$beta_true,
                    sigma_g2 = sg, block_size = ANCHOR$block_size, seed = seed)
  res <- run_methods(d, methods, K_max = ANCHOR$K_max,
                     theta = ANCHOR$theta, N_delta = ANCHOR$N_delta,
                     tau2 = ANCHOR$tau2)
  metrics <- compute_metrics(res, d$truth)
  metrics$n <- ANCHOR$n; metrics$m <- m; metrics$rho <- ANCHOR$rho
  metrics$sigma_g2 <- sg; metrics$rep <- rep_id
  save_cell_rep(out_dir, cell_tag, rep_id,
                  list(metrics = metrics, truth = d$truth, results = res))
  "done"
}

# 06: whole-genome mice — single joint analysis over the full marker panel.
# Implementation lives in `06_realdata_mice_wholegenome.R`.
run_cell_06mice_whole <- function(dat, methods,
                                     K_max = 30L, theta = 0.99,
                                     N_delta = 10L) {
  out_dir <- "results/bench_full/06_realdata_mice_wholegenome"
  out_file <- file.path(out_dir, "whole_genome.rds")
  if (file.exists(out_file) && file.size(out_file) > 0) return("skip")
  d <- list(X = dat$X, y = dat$y, truth = NULL)
  res <- run_methods(d, methods, K_max = K_max, theta = theta,
                     N_delta = N_delta, tau2 = ANCHOR$tau2)

  # Decorate with target-chromosome info
  CHR_LEVELS <- c(as.character(1:19), "X")
  for (m_name in names(res)) {
    sel <- res[[m_name]]$indices
    if (length(sel) > 0L) {
      res[[m_name]]$snp_names  <- dat$map$snp_id[sel]
      res[[m_name]]$snp_chrs   <- as.character(dat$map$chr[sel])
      cnt <- table(factor(res[[m_name]]$snp_chrs, levels = CHR_LEVELS))
      res[[m_name]]$chr_counts <- as.integer(cnt)
      names(res[[m_name]]$chr_counts) <- CHR_LEVELS
      res[[m_name]]$target_chrs <- CHR_LEVELS[cnt > 0]
    } else {
      res[[m_name]]$snp_names  <- character(0)
      res[[m_name]]$snp_chrs   <- character(0)
      res[[m_name]]$chr_counts <- setNames(rep(0L, length(CHR_LEVELS)),
                                              CHR_LEVELS)
      res[[m_name]]$target_chrs <- character(0)
    }
  }

  chr_matrix <- do.call(cbind, lapply(res, function(r) r$chr_counts))
  rownames(chr_matrix) <- CHR_LEVELS
  consensus_count <- rowSums(chr_matrix > 0)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(
    n = nrow(dat$X), m = ncol(dat$X), map = dat$map,
    chr_levels = CHR_LEVELS, chr_matrix = chr_matrix,
    consensus = list(
      hit_count_per_chr = setNames(consensus_count, CHR_LEVELS),
      chr_hit_by_at_least_2 = CHR_LEVELS[consensus_count >= 2],
      chr_hit_by_at_least_3 = CHR_LEVELS[consensus_count >= 3],
      chr_hit_by_at_least_4 = CHR_LEVELS[consensus_count >= 4]),
    results = res), out_file)
  "done"
}

# 05: chromosome-wise mice (data loaded once, passed to mclapply via closure)
run_cell_05mice <- function(chr_id, dat, methods,
                              K_max = 5L, theta = 0.95, N_delta = 10L) {
  out_dir <- "results/bench_full/05_realdata_mice"
  out_file <- file.path(out_dir, sprintf("chr%s.rds", chr_id))
  if (file.exists(out_file) && file.size(out_file) > 0) return("skip")
  cols <- which(dat$map$chr == chr_id)
  if (length(cols) < 5L) return("too_small")
  X_chr  <- dat$X[, cols, drop = FALSE]
  snp_id <- dat$map$snp_id[cols]
  d <- list(X = X_chr, y = dat$y, truth = NULL)
  res <- run_methods(d, methods, K_max = K_max, theta = theta,
                     N_delta = N_delta, tau2 = ANCHOR$tau2)
  for (m_name in names(res)) {
    sel <- res[[m_name]]$indices
    res[[m_name]]$snp_names <- if (length(sel) > 0L) snp_id[sel] else character(0)
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(chr = chr_id, m = ncol(X_chr), n = nrow(X_chr),
                 snp_id = snp_id, results = res), out_file)
  "done"
}


# =============================================================================
# Build full task list (one task = one (axis, cell, rep) atomic unit)
#
# `demo = TRUE` returns a tiny subset of one (cheapest) cell per axis,
# B forced to 2, plus a single small chromosome for mice. Total wall-clock
# on a 4-core laptop should stay under ~5 min and is intended for
# pipeline smoke-testing rather than science.
# =============================================================================
build_tasks <- function(B, axes, demo = FALSE) {
  if (demo) {
    # Demo grid: smallest cell per axis, B=2, only chr 19 (smallest, 249 SNPs)
    # Demo overrides ANCHOR with smaller m and tighter K_max via
    # ANCHOR_demo (set in run_all() when demo=TRUE).
    n_grid     <- 300L                # smaller than ANCHOR$n
    rho_grid   <- 0.95
    sig_grid   <- "medium"
    m_grid     <- 1000L               # 5x smaller than ANCHOR$m
    chr_grid   <- "19"
    sg_grid    <- 0.0
    B          <- min(B, 2L)
  } else {
    n_grid   <- c(500L, 1000L, 3000L)
    rho_grid <- c(0.80, 0.95, 0.98)
    sig_grid <- c("weak", "medium", "strong")
    m_grid   <- c(5000L, 10000L)
    chr_grid <- c(as.character(1:19), "X")
    sg_grid  <- c(0.0, 0.5)
  }

  tasks <- list()

  if ("01_n" %in% axes) {
    for (n in n_grid) for (sg in sg_grid) for (rep_id in seq_len(B))
      tasks[[length(tasks) + 1L]] <- list(
        axis = "01_n", n = n, sg = sg, rep_id = rep_id,
        label = sprintf("01n[n=%d sg=%.1f b%02d]", n, sg, rep_id))
  }

  if ("02_rho" %in% axes) {
    for (rho in rho_grid) for (sg in sg_grid) for (rep_id in seq_len(B))
      tasks[[length(tasks) + 1L]] <- list(
        axis = "02_rho", rho = rho, sg = sg, rep_id = rep_id,
        label = sprintf("02r[rho=%.2f sg=%.1f b%02d]", rho, sg, rep_id))
  }

  if ("03_arch" %in% axes) {
    for (signal in sig_grid) for (sg in sg_grid) for (rep_id in seq_len(B))
      tasks[[length(tasks) + 1L]] <- list(
        axis = "03_arch", signal = signal, sg = sg, rep_id = rep_id,
        label = sprintf("03a[%s sg=%.1f b%02d]", signal, sg, rep_id))
  }

  if ("04_m" %in% axes) {
    for (m in m_grid) for (sg in sg_grid) for (rep_id in seq_len(B))
      tasks[[length(tasks) + 1L]] <- list(
        axis = "04_m", m = m, sg = sg, rep_id = rep_id,
        label = sprintf("04m[m=%d sg=%.1f b%02d]", m, sg, rep_id))
  }

  if ("05_mice" %in% axes) {
    for (chr_id in chr_grid)
      tasks[[length(tasks) + 1L]] <- list(
        axis = "05_mice", chr = chr_id,
        label = sprintf("05mice[chr%s]", chr_id))
  }

  if ("06_mice_whole" %in% axes && !demo) {
    tasks[[length(tasks) + 1L]] <- list(
      axis = "06_mice_whole",
      label = "06mice_whole[all chromosomes pooled]")
  }

  tasks
}



# one function that knows how to run any task

run_one_task <- function(task, methods, mice_dat = NULL) {
  status <- tryCatch({
    switch(task$axis,
      "01_n"   = run_cell_01n(task$n, task$sg, task$rep_id, methods),
      "02_rho" = run_cell_02rho(task$rho, task$sg, task$rep_id, methods),
      "03_arch"= run_cell_03arch(task$signal, task$sg, task$rep_id, methods),
      "04_m"   = run_cell_04m(task$m, task$sg, task$rep_id, methods),
      "05_mice"= run_cell_05mice(task$chr, mice_dat, methods),
      "06_mice_whole" = run_cell_06mice_whole(mice_dat, methods),
      stop("Unknown axis: ", task$axis))
  }, error = function(e) sprintf("ERROR: %s", conditionMessage(e)))
  list(label = task$label, status = status)
}



# Top-level orchestrator

#' Run the full benchmark suite locally on `n_cores` cores.
#'
#' @param n_cores integer; number of parallel workers (mclapply).
#' @param B integer; replicates per simulation cell.
#' @param axes character vector; subset of
#'   c("01_n","02_rho","03_arch","04_m","05_mice","06_mice_whole").
#' @param with_bslmm logical; include BSLMM (slow MCMC) in method set.
#' @param verbose logical; print one line per completed task.
#' @param demo logical; tiny subset (1 cell/axis, B=2, 1 chromosome, no
#'   sigma_g2 sweep). Designed for ~3-5 min smoke testing.
#' @return invisible; produces RDS checkpoints + (optionally) summary CSVs.
run_all <- function(n_cores = max(1L, detectCores() - 2L),
                      B = 30L,
                      axes = c("01_n", "02_rho", "03_arch", "04_m",
                                "05_mice", "06_mice_whole"),
                      with_bslmm = FALSE,
                      verbose = TRUE,
                      aggregate_at_end = TRUE,
                      demo = FALSE) {

  # Make output paths work regardless of CWD: chdir to project root
  # for the duration of the run, restore on exit.
  saved_wd <- getwd()
  setwd(PROJECT_ROOT)
  on.exit(setwd(saved_wd), add = TRUE)

  # In demo mode, default to including BSLMM (small grid, BSLMM stays fast).
  # In full mode, BSLMM stays opt-in (slow MCMC at large n, m).
  if (demo && missing(with_bslmm)) with_bslmm <- TRUE

  methods <- default_methods(include_bslmm = with_bslmm)
  if (demo) {
    cat("\n=== bench_full DEMO MODE ===\n")
    cat("  Tiny grid (1 cell per axis, B=2, chr 19 only, sigma_g2=0)\n")
    cat("  Anchor overrides: n=300, m=1000, K_max=5, N_delta=10\n")
    cat(sprintf("  Methods: %s\n", paste(methods, collapse = ", ")))
    cat("  Expected wall-clock: 5-15 min on 4 cores (BLAS threads = 1)\n")
    cat("  Use run_all(demo=FALSE) for the full benchmark.\n")
    cat("  Use run_all(demo=TRUE, with_bslmm=FALSE) to skip BSLMM.\n")
    # Override ANCHOR with smaller demo values; mclapply forks see this.
    ANCHOR_save <- ANCHOR
    ANCHOR$m       <- 1000L
    ANCHOR$K_max   <- 5L
    ANCHOR$N_delta <- 10L
    assign("ANCHOR", ANCHOR, envir = .GlobalEnv)
    on.exit(assign("ANCHOR", ANCHOR_save, envir = .GlobalEnv), add = TRUE)
  }
  cat(sprintf("\n=== bench_full: %d cores, B=%d, axes=%s%s ===\n",
              n_cores, if (demo) min(B, 2L) else B,
              paste(axes, collapse = ","),
              if (demo) " [DEMO]" else ""))
  cat(sprintf("Methods: %s\n", paste(methods, collapse = ",")))

  # Pre-load mice if needed (workers inherit it via fork copy-on-write)
  mice_dat <- NULL
  if (any(c("05_mice", "06_mice_whole") %in% axes)) {
    cat("Loading mice data ...\n")
    if (!requireNamespace("BGLR", quietly = TRUE))
      stop("BGLR required for 05_mice / 06_mice_whole. install.packages('BGLR').")
    .env <- new.env()
    data(mice, package = "BGLR", envir = .env)
    X_full <- .env$mice.X; pheno <- .env$mice.pheno; map <- .env$mice.map
    y_full <- pheno$Obesity.BMI
    ok <- which(is.finite(y_full))
    X_ok <- X_full[ok, ]; y_ok <- y_full[ok]
    af <- colMeans(X_ok, na.rm = TRUE) / 2; maf <- pmin(af, 1 - af)
    keep <- which(is.finite(maf) & maf > 0.05)
    mice_dat <- list(X = scale(X_ok[, keep]), y = y_ok,
                       map = map[keep, , drop = FALSE])
    cat(sprintf("  After QC: n=%d, m=%d, %d chromosomes\n",
                nrow(mice_dat$X), ncol(mice_dat$X),
                length(unique(mice_dat$map$chr))))
  }

  tasks <- build_tasks(B, axes, demo = demo)
  n_tasks <- length(tasks)
  cat(sprintf("Total tasks: %d (already-done tasks will be skipped)\n\n",
              n_tasks))

  t_start <- Sys.time()
  if (n_cores > 1L && .Platform$OS.type == "unix") {
    # Fork-based parallel (macOS / Linux)
    results <- mclapply(tasks, function(task) {
      out <- run_one_task(task, methods = methods, mice_dat = mice_dat)
      if (verbose) message(sprintf("  [%s] %s", out$status, out$label))
      out
    }, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    if (n_cores > 1L)
      message("Note: mclapply unavailable on Windows; running serially.")
    results <- vector("list", n_tasks)
    for (i in seq_along(tasks)) {
      results[[i]] <- run_one_task(tasks[[i]], methods = methods,
                                      mice_dat = mice_dat)
      if (verbose) cat(sprintf("  [%4d/%d %s] %s\n",
                                  i, n_tasks, results[[i]]$status,
                                  results[[i]]$label))
    }
  }
  t_elapsed <- as.numeric(Sys.time() - t_start, units = "secs")

  # Status summary
  n_done <- sum(vapply(results, function(r) r$status == "done", logical(1)))
  n_skip <- sum(vapply(results, function(r) r$status == "skip", logical(1)))
  n_err  <- sum(vapply(results, function(r)
    grepl("^ERROR", r$status), logical(1)))
  cat(sprintf("\n=== bench_full done in %.0f s (%.1f min) ===\n",
              t_elapsed, t_elapsed / 60))
  cat(sprintf("  done: %d   skipped (already done): %d   errors: %d\n",
              n_done, n_skip, n_err))
  if (n_err > 0L) {
    cat("  Errors:\n")
    for (r in results) if (grepl("^ERROR", r$status))
      cat(sprintf("    %s -> %s\n", r$label, r$status))
  }

  if (aggregate_at_end) {
    cat("\nAggregating ...\n")
    # In demo mode, single-cell-per-axis means figures degenerate to a
    # single point per facet; that's OK -- plot_axis handles it gracefully.
    agg_path <- file.path(BENCH_DIR, "99_aggregate.R")
    source(agg_path, local = TRUE)
  }
  invisible(results)
}



# CLI entry point

if (sys.nframe() == 0L) {
  args <- parse_args()
  n_cores    <- if (!is.null(args[["n-cores"]])) as.integer(args[["n-cores"]])
                else max(1L, parallel::detectCores() - 2L)
  B          <- if (!is.null(args$B)) as.integer(args$B) else 30L
  with_bslmm <- isTRUE(args[["with-bslmm"]])
  demo       <- isTRUE(args[["demo"]])
  axes_arg   <- if (!is.null(args$axes)) args$axes else "all"
  axes <- if (axes_arg == "all")
             c("01_n", "02_rho", "03_arch", "04_m",
               "05_mice", "06_mice_whole")
          else strsplit(axes_arg, ",")[[1]]

  run_all(n_cores = n_cores, B = B, axes = axes,
          with_bslmm = with_bslmm, demo = demo)
}
