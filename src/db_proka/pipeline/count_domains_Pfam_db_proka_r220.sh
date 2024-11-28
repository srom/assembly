#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -e error_run_domain_counts_Pfam_qsub.txt
#PBS -o output_run_domain_counts_Pfam_qsub.txt

set -e

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220
HMM_DB=/rds/general/project/lms-warnecke-raw/live/Pfam-A_37/Pfam-A.hmm

# Input variables can be set with qsub -v
# If not set the following defaults are used:
if [ -z "${BASE_FOLDER}" ]; then
    BASE_FOLDER=$DB_PROKA_BASE
fi
if [ -z "${METADATA}" ]; then
    METADATA="${DB_PROKA_BASE}/db_proka_r220_genus_10s_metadata.csv"
fi
if [ -z "${OUTPUT}" ]; then
    OUTPUT="${DB_PROKA_BASE}/db_proka_r220_genus_10s_Pfam-A_summary.tsv"
fi

cd $HOME
. load_conda.sh

cd ${HOME}/assembly-main

python -m src.postprocessing.count_domains \
	-i $BASE_FOLDER \
	-d $HMM_DB \
    -o $OUTPUT \
    --metadata_path $METADATA \
	--n_cpus 64 \
    > output_run_domain_counts_Pfam.txt \
    2> error_run_domain_counts_Pfam.txt
