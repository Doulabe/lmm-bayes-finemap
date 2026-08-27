# CBF-LMM: Conditional Bayes-Factor Selection with an Adaptive Polygenic Background

**CBF-LMM** is a conditional Bayes-factor stepwise procedure for **refining**
candidate regions under an **adaptive polygenic background**: previously
selected SNPs enter as fixed effects while the remaining markers define the
linear mixed-model background kernel, rebuilt along the selection path.  At
each step a candidate SNP enters as one additional fixed effect; the residual variance is integrated out analytically, the
variance-component ratio $\delta$ is **numerically marginalised** by
one-dimensional Gauss–Legendre quadrature, and selection is controlled by the
extended BIC.  A fast **plug-in REML** evaluation of $\delta$ and a
posterior-score stopping rule are available as **sensitivity options**, and an
exploratory **joint Schur-complement** score is retained for future development.

> Reference: K. N. Doulabe and L. Lakhal-Chaieb, *Conditional Bayes-Factor
> Selection with an Adaptive Polygenic Background*, 2026 (submitted).

The method is positioned as a **second-stage refinement tool** downstream of a
genome-wide screen or a marginal Bayesian fine-mapper (e.g. SuSiE), not as a
replacement for them: it returns a compact, parsimony-controlled set of
candidates.

## Paper-to-code name mapping

| Paper | Code | Status |
|---|---|---|
| **CBF-LMM** (primary) | `CBF_LMM_stepwise()` = `MS_L_LMM_stepwise_fast()` with defaults (`delta_eval = "marginal"`, `criterion = "eBIC"`) | primary |
| Plug-in REML evaluation | `delta_eval = "reml"` | sensitivity |
| Posterior-score stopping | `criterion = "JointPosterior"` | sensitivity |
| Joint Schur-complement score | `JS_L_LMM_stepwise_fast()` | exploratory / future work |

## Highlights

- **Closed-form conditional Bayes factor** at fixed $\delta$, valid for any
  step $k \ge 1$ under a point-normal $g$-prior on the candidate effect, a
  Jeffreys prior on $\sigma^2$, and weakly-informative priors on the nuisance
  parameters.  Strictly generalises the Wakefield approximate Bayes factor
  (recovered as $\delta \to 0$ at $k = 1$).
- **Primary operating point (CBF-LMM)**: $\delta$-marginalised evaluation
  (Gauss–Legendre quadrature over a half-Cauchy prior on $\sqrt\delta$) with
  eBIC stopping ($\gamma = 1$, parsimony-oriented).
- **Sensitivity options**: a fast plug-in **REML** evaluation of $\delta$
  (justified by an $O_p(n^{-1/2})$ asymptotic-equivalence result and in close
  empirical agreement at fine-mapping sample sizes) and a posterior-score
  stopping rule (`criterion = "JointPosterior"`, recall-oriented).
- **Covariate extension**: the conditioning block extends to arbitrary fixed
  covariates (sex, batch, family, …) without modification to the Bayes-factor
  algebra; this extension and the exploratory **joint Schur-complement** score
  (more stable under within-block LD saturation) are discussed as future
  directions in the manuscript.
- **Comparators** (same $X$, $y$, GRM): BSLMM and BayesR (joint Bayesian LMM),
  SuSiE (credible-set), and FaST-LMM (frequentist single-SNP).
- **Sensitivity suite**: slab variance $\tau^2$, a closed-form empirical-Bayes
  slab (EB-$\tau^2$), the eBIC penalty $\gamma$, the prior on $\delta$
  (half-Cauchy vs inverse-gamma), and a $\delta$-profile identifiability
  diagnostic on the mouse panel.
- **Reproducibility**: full benchmark pipeline with seed-deterministic
  cell-level checkpointing.

## Quick start

```r
# Load the framework
source("R/LMM_core.R")
source("R/LMM_stepwise_fast.R")
source("R/LMM_reml_bf.R")        # needed for the delta_eval = "reml" sensitivity option
source("R/CBF_LMM.R")            # paper-facing alias

# Simulate a small dataset: n=500, m=300, K_true=3 active SNPs
n <- 500L; m <- 300L; K_true <- 3L
set.seed(42)
X <- matrix(rnorm(n * m), n, m)
X <- scale(X)
truth <- c(50L, 120L, 200L)
beta_true <- c(0.8, 0.4, 0.4)
y <- as.numeric(X[, truth] %*% beta_true) + rnorm(n)

# Primary procedure (CBF-LMM): marginalised delta + eBIC stopping — the defaults.
res <- CBF_LMM_stepwise(y, X, tau2 = 0.04, K_max = 10L, n_nodes = 15L)
cat("Selected indices:", res$indices, "\n")
cat("Truth:           ", truth, "\n")
cat("K_hat:", res$K_hat, " (true K_true =", K_true, ")\n")

# Sensitivity option: plug-in REML evaluation of delta
res_reml <- CBF_LMM_stepwise(y, X, tau2 = 0.04, K_max = 10L,
                             delta_eval = "reml")
```

See [`examples/quickstart.R`](examples/quickstart.R) for a complete runnable
example and [`examples/real_data_example.R`](examples/real_data_example.R) for a
minimal mouse-BMI analysis.

### Choosing the slab variance $\tau^2$ (empirical-Bayes option)

The default $\tau^2 = 0.04$ corresponds to a prior standard deviation of $0.2$
residual-scale units for a standardized candidate effect; in applications it
can be anchored to the anticipated standardized effect size.  Because the
conditional Bayes factor is available in closed form in $\tau^2$, a per-step
**empirical-Bayes estimate** of the slab along the selection path is provided:

```bash
Rscript sim/bench_full/21_eb_tau2_sensitivity.R
```

The script ranks candidates at the default, maximizes the closed-form Bayes
factor in $\tau^2$ for the leading candidate at each step (`tau2_mode = "eb"`
in the in-script stepwise), reports the distribution of the per-step estimates,
and verifies that the selected sets match the fixed default on the anchor
design (Supplementary Note S5 of the manuscript).

## Repository layout

```
.
├── README.md                  This file
├── INSTALL.md                 Detailed installation / dependency notes
├── LICENSE                    MIT
├── CITATION.cff               Citation metadata
├── R/                         Core framework implementations
│   ├── LMM_core.R                Closed-form Bayes-factor primitives
│   ├── LMM_stepwise.R            Per-candidate reference implementation
│   ├── LMM_stepwise_fast.R       Batched version (delta_eval = "marginal"/"reml")
│   ├── LMM_reml_bf.R             Plug-in REML-BF evaluation (delta_eval = "reml")
│   ├── LMM_stepwise_mixtau.R     Mixture-slab variant (sensitivity)
│   ├── LMM_stepwise_pathavg.R    Path-averaged variant (exploratory)
│   └── CBF_LMM.R                 Paper-facing alias CBF_LMM_stepwise()
├── sim/                       Simulation and benchmark scripts
│   ├── bench_full/               Full reproducibility pipeline
│   │   ├── 00_config.R              Common configuration + comparator runners
│   │   ├── 01-04_*.R                Four-axis benchmark (n, ρ, signal, m)
│   │   ├── 01b-04b_*_remlbf.R       Same four axes under the REML-BF evaluation
│   │   ├── 05-06_*.R                Mouse BMI (chromosome-wise + whole-genome)
│   │   ├── 07_prior_sensitivity.R   Half-Cauchy vs inverse-gamma on δ
│   │   ├── 08_add_fastlmm.R         Add FaST-LMM comparator to existing RDS
│   │   ├── 10_threshold_independent.R, 10c_*  Threshold-independent ranking
│   │   ├── 11_tau2_sensitivity.R    τ² sensitivity
│   │   ├── 12_mouse_autosomes_*.R   Autosomes-only mouse BMI scan
│   │   ├── 13_mixtau_benchmark.R    Mixture-slab benchmark
│   │   ├── 14_pathavg_benchmark.R   Path-averaging benchmark (exploratory)
│   │   ├── 15_mouse_bmi_with_sex.R  Mouse BMI with sex covariate
│   │   ├── 15b_mouse_reml_bf_sensitivity.R  REML-BF real-data sensitivity scan
│   │   ├── 16_ablation.R            Variance-ratio × scoring × stopping ablation
│   │   ├── 17_rho_V_diagnostic.R    Empirical ρ_V / D_jk distribution under LD
│   │   ├── 18_bootstrap_CI.R        Bootstrap confidence intervals
│   │   ├── 19_threshold_full_grid.R Top-K* recall on the full 16-cell grid
│   │   ├── 21_eb_tau2_sensitivity.R Closed-form empirical-Bayes slab (EB-τ²)
│   │   ├── 22_gamma_sensitivity.R   eBIC penalty γ sensitivity
│   │   ├── 23_delta_profile_mouse.R δ-profile identifiability diagnostic
│   │   ├── 25_semisynth_1000g.R    Semi-synthetic benchmark on real 1000G panels
│   │   ├── 99_aggregate.R           Pool & summarise raw RDS files
│   │   ├── make_*.R                 LaTeX table / figure generators
│   │   ├── make_cbf_figures.R       CBF-LMM manuscript figures (5-method roster)
│   │   ├── make_cbf_rho_table.R     CBF-LMM LD-axis supplement table
│   │   ├── 26_geuvadis_eqtl.R       GEUVADIS human cis-eQTL illustration
│   │   ├── make_geuvadis_table.R    GEUVADIS results table
│   │   └── run_all.R                Single-command driver (axes 01–06)
│   └── validate_lmm_bayes.R      Six-check internal validation (V1–V6)
├── examples/                  Minimal demos
│   ├── quickstart.R              20-line single-dataset demo
│   ├── simulation_example.R      One anchor cell with full pipeline
│   ├── real_data_example.R       Minimal mouse-BMI analysis
│   └── delta_profile_diagnostic.R  δ-profile QC plot (REML-BF vs MBF rule)
└── tests/                     Sanity checks (R)
```

## Reproducing the benchmark

The simulation grid has **16 unique cells** spanning four axes — sample size
$n \in \{500, 1000, 3000\}$, block-AR(1) correlation $\rho \in \{0.80, 0.95,
0.98\}$, signal architecture $\in \{$weak, medium, strong$\}$, and number of
SNPs $m \in \{5000, 10000\}$ — each crossed with the homogeneous-noise and
polygenic regimes $\sigma_g^2 \in \{0, 0.5\}$.  Everything reproduces from
([`sim/bench_full/`](sim/bench_full/)):

```bash
# 1.  Core four-axis grid (MBF evaluation) + aggregation
Rscript sim/bench_full/run_all.R
Rscript sim/bench_full/99_aggregate.R --filter-n
Rscript sim/bench_full/make_aggregate_table.R
Rscript sim/bench_full/make_runtime_table.R
Rscript sim/bench_full/make_pub_figures.R
Rscript sim/bench_full/make_supp_tables.R

# 2.  REML-BF evaluation on the same four axes (lead row of the aggregate table)
Rscript sim/bench_full/01b_scaling_n_remlbf.R
Rscript sim/bench_full/02b_scaling_rho_remlbf.R
Rscript sim/bench_full/03b_arch_remlbf.R
Rscript sim/bench_full/04b_scaling_m_remlbf.R

# 3.  Ablation (variance-ratio × scoring × stopping) and bootstrap intervals
Rscript sim/bench_full/16_ablation.R
Rscript sim/bench_full/18_bootstrap_CI.R

# 4.  Threshold-independent ranking: Top-K* recall on the full 16-cell grid
Rscript sim/bench_full/10_threshold_independent.R
Rscript sim/bench_full/19_threshold_full_grid.R
Rscript sim/bench_full/make_threshold_indep_table.R
Rscript sim/bench_full/make_threshold_strat_tables.R   # by-axis stratified tables

# 5.  Prior / hyperparameter sensitivity
Rscript sim/bench_full/07_prior_sensitivity.R      # half-Cauchy vs inverse-gamma on δ
Rscript sim/bench_full/make_sensitivity_table.R
Rscript sim/bench_full/11_tau2_sensitivity.R       # slab variance τ²
Rscript sim/bench_full/make_tau2_table.R
Rscript sim/bench_full/21_eb_tau2_sensitivity.R    # empirical-Bayes slab EB-τ²
Rscript sim/bench_full/22_gamma_sensitivity.R      # eBIC penalty γ

# 6.  Mouse BMI: whole-genome ± sex × autosomes-only, REML-BF scan, δ-profile
Rscript sim/bench_full/06_realdata_mice_wholegenome.R
Rscript sim/bench_full/12_mouse_autosomes_only.R
Rscript sim/bench_full/15_mouse_bmi_with_sex.R
Rscript sim/bench_full/15b_mouse_reml_bf_sensitivity.R
Rscript sim/bench_full/23_delta_profile_mouse.R
Rscript sim/bench_full/make_autosomes_table.R
Rscript sim/bench_full/make_bmi_sex_table.R

# 7.  Semi-synthetic benchmark on real 1000G panels (main text Table 5)
#     Requires the two locus panels — see "Semi-synthetic 1000G benchmark" below.
Rscript sim/bench_full/25_semisynth_1000g.R --B 200
Rscript sim/bench_full/make_semisynth_tables.R         # pooled + per-locus tables

# 8.  CBF-LMM manuscript figures and LD-axis supplement table
Rscript sim/bench_full/make_cbf_figures.R
Rscript sim/bench_full/make_cbf_rho_table.R

# 9.  GEUVADIS human cis-eQTL illustration (see data section below)
Rscript sim/bench_full/26_geuvadis_eqtl.R
Rscript sim/bench_full/make_geuvadis_table.R
```

### Semi-synthetic 1000G benchmark

Step 7 evaluates the same eight methods on **real** genotypes (real LD) rather
than block-AR(1) simulations.  The two locus panels are **not redistributed
here**; they are built from the public 1000 Genomes phase-3 release.  Each panel
is an `.rds` holding a list with

| field | content |
|---|---|
| `G` | `n × m` genotype matrix, integer dosages `0/1/2` |
| `R` | `m × m` LD matrix (`cor(G)`) |
| `pos` | physical positions, length `m` |
| `snpvar` | per-SNP prior variance (unused by this script) |

The published run uses the EUR subset (`n = 503`, MAF ≥ 0.05) at a chromosome-1
locus (`m = 1,493`) and a chromosome-6 locus (`m = 2,143`).  Place them at
`data/1000g/locus_chr1_1000g.rds` and `data/1000g/locus_chr6_1000g.rds`, or point
the script elsewhere:

```bash
Rscript sim/bench_full/25_semisynth_1000g.R \
    --locus_chr1 /path/to/chr1.rds --locus_chr6 /path/to/chr6.rds --B 200
```

Recovery is scored in the LD-aware sense standard for real-LD fine-mapping (a
causal counts as recovered, and a selection as a true positive, at
`r² ≥ 0.25`), because at `n = 503` a causal variant routinely has a perfect-LD
twin and exact-index recovery is not identifiable.

### GEUVADIS cis-eQTL illustration

Step 9 fine-maps eight strong cis-eQTL genes in the GEUVADIS LCL panel
(n = 358 individuals shared with the 1000 Genomes phase-3 EUR reference),
using the published EUR373 best-association list as a concordance benchmark.
All inputs are public; genotypes are streamed remotely per locus by
`bcftools` (required on PATH):

```bash
mkdir -p data/geuvadis
B=http://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/analysis_results
curl -o data/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz \
     $B/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
curl -o data/geuvadis/EUR373.gene.cis.FDR5.best.rs137.txt.gz \
     $B/EUR373.gene.cis.FDR5.best.rs137.txt.gz
# sample overlap: GEUVADIS columns intersected with the phase-3 EUR panel
curl -s http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel \
  | awk '$3=="EUR"{print $1}' | sort > data/geuvadis/eur_phase3.txt
gzcat data/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz | head -1 \
  | tr '\t' '\n' | grep -E '^HG|^NA' | sort > data/geuvadis/geuvadis_samples.txt
comm -12 data/geuvadis/eur_phase3.txt data/geuvadis/geuvadis_samples.txt \
  > data/geuvadis/geuvadis_eur_overlap.txt   # 358 individuals
```

Headline result (best selected variant; **bold** tags the published lead at
r2 = 1.00; "---" = empty selection at the method's prespecified threshold):

| Gene | Lead eQTL | CBF-LMM | SuSiE | BSLMM | BayesR | FaST-LMM |
|---|---|---|---|---|---|---|
| ERAP2 | rs2910686 | **rs2927608** | --- | **rs2910686** | **rs2910686** | --- |
| RPS26 | rs10876864 | **rs10876864** | --- | **rs10876864** | **rs10876864** | --- |
| SLFN5 | rs11080327 | **rs11080327** | **rs11080327** | rs883416 (0.97) | **rs11080327** | **rs11080327** |
| SNHG5 | rs1059307 | **rs1059307** | **rs1059307** | **rs1059307** | **rs1059307** | --- |
| FLVCR1-AS1 | rs12123978 | **rs61832055** | --- | **rs11120042** | **rs10864005** | --- |
| PEX6-region | rs6907751 | rs9986447 +1 (0.88) | **rs6907751** +1 | rs2296804 (0.63) | rs2296805 (0.63) | --- |
| TRA2A-AS | rs10233039 | **rs6461691** | --- | rs10266123 (0.99) | --- | --- |
| ZNF266 | rs10420709 | **rs11878970** | --- | **rs10411141** | --- | --- |

CBF-LMM returns a selection at all eight loci and tags the published lead at
r2 = 1.00 at seven of eight; SuSiE's PIP > 0.99 selection is empty at five
loci (posterior mass split across perfect-LD proxies).


Each script supports `--cores N` for `mclapply` parallelism and writes per-cell
RDS checkpoints; re-running an interrupted job resumes from the last completed
checkpoint.  **Important:** set
`OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1` to avoid Accelerate/OpenBLAS thread contention under
`mclapply` workers (the most common source of order-of-magnitude slowdowns).

> Note: the comparison is restricted to methods that model the polygenic
> background jointly (BSLMM, BayesR) or via credible-set decomposition
> (SuSiE), plus a frequentist LMM reference (FaST-LMM).  A marginal
> non-LMM spike-and-slab (`varbvs`) falls outside this comparator class
> and is not included.

## Internal validation

A six-check validation suite is provided in
[`sim/validate_lmm_bayes.R`](sim/validate_lmm_bayes.R).  V1 confirms
machine-precision recovery of Wakefield's ABF at $\delta \to 0$; V2 checks
asymptotic agreement with plug-in REML at $n = 1{,}000$; V3–V6 cover type-I
behaviour under the null, quadrature convergence, stepwise recall, and the
joint-vs-marginal distinction.

```bash
Rscript sim/validate_lmm_bayes.R
```

## Installation

See [`INSTALL.md`](INSTALL.md).  Briefly:

```r
install.packages(c("statmod", "Matrix",                 # core
                   "dplyr", "tidyr", "ggplot2",         # tables/figures
                   "susieR", "rrBLUP", "BGLR",          # comparators + mouse data
                   "hibayes"))                          # BayesR comparator
```

All packages are on CRAN; no GitHub-only dependencies.

## Citation

If you use this code in your work, please cite the manuscript (see
[`CITATION.cff`](CITATION.cff)):

```
@article{Doulabe2026LMMBayes,
  author  = {Doulabe, Kossi N. and Lakhal-Chaieb, Lajmi},
  title   = {Conditional Bayes-Factor Selection with an Adaptive
             Polygenic Background},
  journal = {(submitted)},
  year    = {2026}
}
```

## License

MIT — see [`LICENSE`](LICENSE).
