# Results directory

Run the benchmark scripts in `sim/bench_full/` to populate this directory.

Each script writes per-(cell, replicate) RDS files to subdirectories
named after the script ID (e.g., `01_scaling_n/`, `13_mixtau/`).

The aggregation scripts (`99_aggregate.R`, `make_*.R`) read these RDS
files and write summary CSVs and LaTeX tables.

This directory is gitignored except for this README.

