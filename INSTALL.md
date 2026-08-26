# Installation

## System requirements

- **R** ≥ 4.0  (tested on R 4.4 arm64 / x86_64)
- A C/C++ compiler (for the `Rcpp` dependency of some comparator
  packages); on macOS install Xcode command-line tools, on Linux any
  recent gcc/clang.
- Enough RAM for an `n × n` GRM eigendecomposition: ~8 GB suffices for
  $n \le 3{,}000$ used in the benchmark.  For the mouse BMI application
  ($n = 1{,}814$, $m = 10{,}346$): ~16 GB recommended.

## R package dependencies

### Core (always required)

```r
install.packages(c("statmod", "Matrix"))
```

- `statmod` — Gauss–Legendre quadrature nodes/weights.
- `Matrix` — sparse-matrix utilities used by the eigendecomposition path.

### Tables and figures

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "scales"))
```

### External comparator methods (benchmark)

```r
install.packages(c("susieR", "rrBLUP", "BGLR", "hibayes"))
```

- `susieR` — Sum-of-Single-Effects credible-set comparator (Wang et al.,
  2020).
- `rrBLUP` — REML estimator used by the FaST-LMM-equivalent comparator
  and by the plug-in REML sensitivity evaluation of delta.
- `BGLR` — provides the BSLMM comparator (`model = "BSLMM"`) and the
  Bayesian-Lasso runner; also ships the Valdar et al. (2006)
  heterogeneous-stock mouse panel as `data(mice)`.
- `hibayes` — provides the BayesR comparator (`ibrm(method = "BayesR")`).

All packages are on CRAN; no GitHub-only dependencies.

## BLAS thread control (critical on macOS)

The framework uses `parallel::mclapply` for cross-replicate
parallelism with one thread per worker.  Accelerate (macOS) and
OpenBLAS (Linux) default to multi-threaded BLAS, which **collides
catastrophically** with `mclapply` and produces ~10–100× slowdowns
plus silent worker crashes.  **Always pin to one thread per worker
before launching any benchmark script**:

```bash
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

These environment variables must be set in the shell that invokes
`Rscript`, not in the script itself.

## Quick verification

After installing dependencies, run the validation suite to confirm
the framework is operational:

```bash
Rscript sim/validate_lmm_bayes.R
```

Expected output: six checks (V1–V6) all reporting `PASS`.  Total
runtime $\sim 5$ minutes on a single CPU core.

## Reproducing the manuscript benchmark

See [`README.md`](README.md) section *Reproducing the benchmark*.
Total compute: $\sim 400$ CPU-hours for the $B = 100$ replicates on
$18$ axis cells, executed in $\sim 2$ wall-clock days on a 20-core
host with the BLAS thread settings above.

For a quick demo (single anchor cell, $B = 10$ replicates), use:

```bash
Rscript sim/bench_full/01_scaling_n.R --n 1000 --sg 0 --B 10
```

Output goes to `results/bench_full/01_scaling_n/`.

## Common issues

| Symptom | Cause | Fix |
|---|---|---|
| `Rscript` jobs hang at ~0% CPU | Accelerate BLAS thread contention | Set `OMP_NUM_THREADS=1` etc. (see above) |
| Silent failures, `K_hat = 0` everywhere | `requireNamespace("rrBLUP")` race in workers | Pre-load `rrBLUP` in parent process before `mclapply` (already done in `08b_*.R` and `15_*.R`) |
| `eigen()` errors at large $n$ | Insufficient RAM | Reduce `K_max` or use a machine with more RAM |
| `parse error` in `label_parsed` | Old ggplot2 (<3.4) | Update: `install.packages("ggplot2")` |
| Missing `BGLR::data(mice)` | BGLR version too old | Install latest CRAN: `install.packages("BGLR")` |
