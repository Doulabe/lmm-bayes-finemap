# Tests

Sanity-check and validation scripts for the framework.

## Internal validation suite

[`../sim/validate_lmm_bayes.R`](../sim/validate_lmm_bayes.R) runs
six checks (V1–V6, see manuscript Appendix A):

| Check | Verifies |
|---|---|
| V1 | Closed-form Bayes factor recovers Wakefield's ABF at δ → 0 (to machine precision) |
| V2 | Asymptotic equivalence with plug-in REML at n = 1,000 (Proposition 4.1) |
| V3 | Type-I behaviour under the global null at small n |
| V4 | Quadrature convergence in N_δ |
| V5 | Stepwise recall and stopping at K* = 3 |
| V6 | Joint vs marginal variant difference under block-AR(1) LD |

```bash
Rscript sim/validate_lmm_bayes.R
```

Total runtime: ~5 minutes on a single CPU core.  All six checks should
report `PASS`.

## Mixture-slab sanity check

[`../examples/quickstart.R`](../examples/quickstart.R) and the
mixture-slab wrapper test inside
[`../sim/bench_full/13_mixtau_benchmark.R`](../sim/bench_full/13_mixtau_benchmark.R)
provide a bit-exact identity check: running the mixture wrapper with
a single τ² = 0.04 (weight = 1) must produce the same selected
indices and log-Bayes-factor values as the standard single-τ²
function.

```r
source("R/LMM_stepwise_fast.R")
source("R/LMM_stepwise_mixtau.R")

# ... generate data d ...

res_std <- MS_L_LMM_stepwise_fast(d$y, d$X, tau2 = 0.04, ...)
res_mix <- MS_L_LMM_stepwise_fast_mixtau(d$y, d$X,
                                              tau2_mix = c(0.04),
                                              weights = c(1), ...)

stopifnot(identical(sort(res_std$indices), sort(res_mix$indices)))
stopifnot(max(abs(res_std$log_BF_at_step - res_mix$log_BF_at_step)) < 1e-10)
```

## Notes

The benchmark scripts in `sim/bench_full/` are not formal unit tests
but include explicit consistency checks (e.g., seed-based regeneration
of data matrices to verify framework invariance under restart).  See
the per-script comments for details.
