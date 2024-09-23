#!/bin/bash
set -e

python -m src.fetch_assemblies \
    -l test_data/assembly_accessions.txt \
    -o test_data \
    --cpu 2
