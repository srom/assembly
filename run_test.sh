#!/bin/bash
set -e

python -m src.fetch_assemblies \
    -l test_data/assembly_accessions.txt \
    -o test_data

python -m src.postprocessing.concatenate_proteins \
    -i test_data \
    -o test_data/all_proteins.fasta
