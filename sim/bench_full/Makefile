# =============================================================================
# Makefile -- bench_full orchestrator with GNU parallel.
#
# Run from the project root: `make -f sim/bench_full/Makefile <target>`
# Or:  `cd sim/bench_full && ln -s Makefile ../../Makefile.bench` then
#      `make -f Makefile.bench <target>`.
#
# Targets:
#   all              run scaling_n, scaling_rho, arch, scaling_m, mice (no BSLMM)
#   all-with-bslmm   same + BSLMM (much slower; expect +3h)
#   sims             scaling_n + scaling_rho + arch + scaling_m
#   scaling_n        run script 01 in parallel
#   scaling_rho      run script 02 in parallel
#   arch             run script 03 in parallel
#   scaling_m        run script 04 in parallel
#   mice             run script 05 chromosome-wise in parallel
#   aggregate        collect CSVs and plot figures
#   clean            remove all results/bench_full/ outputs
#   clean-cache      remove only summary CSVs and figures (keep cells)
# =============================================================================

# ---- Knobs ------------------------------------------------------------------
JOBS    ?= 4                  # parallel workers
B       ?= 10                 # replicates per cell
RSCRIPT ?= Rscript
RESULTS  = results/bench_full

# Add `--include-bslmm` here to enable BSLMM in all sim scripts:
BSLMM_FLAG ?=

REP_SEQ := $(shell seq 1 $(B))

# ---- Phony ------------------------------------------------------------------
.PHONY: all all-with-bslmm sims scaling_n scaling_rho arch scaling_m mice \
        aggregate clean clean-cache help

help:
	@echo "Targets:"
	@echo "  all              everything except BSLMM"
	@echo "  all-with-bslmm   same + BSLMM (slow)"
	@echo "  sims             4 simulation axes"
	@echo "  scaling_n        n axis only"
	@echo "  scaling_rho      rho axis only"
	@echo "  arch             signal-strength axis only"
	@echo "  scaling_m        m axis only"
	@echo "  mice             chromosome-wise mouse BMI"
	@echo "  aggregate        collect CSVs + figures"
	@echo "  clean            remove all bench_full outputs"
	@echo ""
	@echo "Variables:"
	@echo "  JOBS=$(JOBS)         parallel workers"
	@echo "  B=$(B)            replicates per cell"
	@echo "  BSLMM_FLAG='$(BSLMM_FLAG)'  set to '--include-bslmm' to enable BSLMM"

# ---- All --------------------------------------------------------------------
all: sims mice aggregate

all-with-bslmm:
	$(MAKE) -f sim/bench_full/Makefile sims mice aggregate \
	  BSLMM_FLAG=--include-bslmm

sims: scaling_n scaling_rho arch scaling_m

# ---- 01: scaling n ----------------------------------------------------------
scaling_n:
	@echo "=== 01_scaling_n: 6 cells x $(B) reps on $(JOBS) workers ==="
	parallel -j $(JOBS) --bar \
	  "$(RSCRIPT) sim/bench_full/01_scaling_n.R --n {1} --sg {2} --rep {3} $(BSLMM_FLAG)" \
	  ::: 500 1000 3000 ::: 0 0.5 ::: $(REP_SEQ)

# ---- 02: scaling rho --------------------------------------------------------
scaling_rho:
	@echo "=== 02_scaling_rho: 6 cells x $(B) reps on $(JOBS) workers ==="
	parallel -j $(JOBS) --bar \
	  "$(RSCRIPT) sim/bench_full/02_scaling_rho.R --rho {1} --sg {2} --rep {3} $(BSLMM_FLAG)" \
	  ::: 0.80 0.95 0.98 ::: 0 0.5 ::: $(REP_SEQ)

# ---- 03: arch ---------------------------------------------------------------
arch:
	@echo "=== 03_arch: 6 cells x $(B) reps on $(JOBS) workers ==="
	parallel -j $(JOBS) --bar \
	  "$(RSCRIPT) sim/bench_full/03_arch.R --signal {1} --sg {2} --rep {3} $(BSLMM_FLAG)" \
	  ::: weak medium strong ::: 0 0.5 ::: $(REP_SEQ)

# ---- 04: scaling m ----------------------------------------------------------
scaling_m:
	@echo "=== 04_scaling_m: 4 cells x $(B) reps on $(JOBS) workers ==="
	parallel -j $(JOBS) --bar \
	  "$(RSCRIPT) sim/bench_full/04_scaling_m.R --m {1} --sg {2} --rep {3} $(BSLMM_FLAG)" \
	  ::: 5000 10000 ::: 0 0.5 ::: $(REP_SEQ)

# ---- 05: mice chromosomes ---------------------------------------------------
mice:
	@echo "=== 05_realdata_mice_full: 20 chromosomes on $(JOBS) workers ==="
	parallel -j $(JOBS) --bar \
	  "$(RSCRIPT) sim/bench_full/05_realdata_mice_full.R --chr {1} $(BSLMM_FLAG)" \
	  ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X

# ---- Aggregate --------------------------------------------------------------
aggregate:
	@echo "=== 99_aggregate ==="
	$(RSCRIPT) sim/bench_full/99_aggregate.R

# ---- Clean ------------------------------------------------------------------
clean:
	rm -rf $(RESULTS)/01_scaling_n $(RESULTS)/02_scaling_rho \
	       $(RESULTS)/03_arch $(RESULTS)/04_scaling_m \
	       $(RESULTS)/05_realdata_mice
	rm -f  $(RESULTS)/*_summary.csv $(RESULTS)/*_raw.csv \
	       $(RESULTS)/fig_*.pdf

clean-cache:
	rm -f $(RESULTS)/*_summary.csv $(RESULTS)/*_raw.csv \
	      $(RESULTS)/fig_*.pdf
