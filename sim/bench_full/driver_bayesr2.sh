#!/bin/bash
# Relaunch of the BayesR long-chain sequence with exit-code checks.
# A failed step prints BAYESR-CHAIN-FAILED and aborts; the success
# sentinel BAYESR-DONE-ALL is printed only if all three steps exit 0.
cd /Users/kossi/Desktop/Dossier_these/Redaction/LMM_Gaussian_Bayesian
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

run_step() {
  local name="$1"; shift
  echo "=== $name $(date) ==="
  "$@"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "BAYESR-CHAIN-FAILED at $name rc=$rc $(date)"
    exit 1
  fi
}

run_step "50 controlled" Rscript sim/bench_full/50_bayesr_long_controlled.R --cores 5
run_step "51 secondary"  Rscript sim/bench_full/51_bayesr_long_secondary.R --cores 5
run_step "52 realdata"   Rscript sim/bench_full/52_bayesr_long_realdata.R
echo "BAYESR-DONE-ALL $(date)"
