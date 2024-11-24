#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=64gb
#PBS -e error_download_db_proka_r220_qsub.txt
#PBS -o output_download_db_proka_r220_qsub.txt

DB_PROKA=/rds/general/project/lms-warnecke-raw/live/db_proka_r220

cd $HOME
. load_conda.sh

cd $HOME/assembly-main

python -m src.fetch_assemblies \
	-l $DB_PROKA/db_proka_r220_genus_10s_accessions.txt \
	-s $DB_PROKA/ncbi_assembly_summary.txt.gz \
	-o $DB_PROKA \
	--cpu=64 \
	> $HOME/output_download_db_proka_r220.txt \
	2> $HOME/error_download_db_proka_r220.txt
