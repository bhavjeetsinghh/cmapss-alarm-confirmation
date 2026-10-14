#!/usr/bin/env bash
# Reproduces every number, table and figure of the paper (about 1 minute).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results figures
(cd src && python cmapss.py ../data/CMAPSSData ../results)
(cd src && python figs4.py ../data/CMAPSSData ../results/cmapss.json ../data/bearing_paper3_results.json ../figures)
echo "Done: results/cmapss.json and figures/"
