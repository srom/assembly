#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -e error_download_db_proka_r220_qsub.txt
#PBS -o output_download_db_proka_r220_qsub.txt

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220

# Input variables can be set with qsub -v
# If not set the following defaults are used:
if [ -z "${ACCESSIONS}" ]; then
    ACCESSIONS="${DB_PROKA_BASE}/db_proka_r220_genus_10s_accessions.txt"
fi
if [ -z "${NCBI_SUMMARY}" ]; then
    NCBI_SUMMARY="${DB_PROKA_BASE}/ncbi_assembly_summary.txt.gz"
fi
if [ -z "${OUTPUT_FOLDER}" ]; then
    OUTPUT_FOLDER=$DB_PROKA_BASE
fi

cd $HOME
. load_conda.sh

cd $HOME/assembly-main

python -m src.fetch_assemblies \
	-l $ACCESSIONS \
	-s $NCBI_SUMMARY \
	-o $OUTPUT_FOLDER \
	--cpu=64 \
	> $HOME/output_download_db_proka_r220.txt \
	2> $HOME/error_download_db_proka_r220.txt
