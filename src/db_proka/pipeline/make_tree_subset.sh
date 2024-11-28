#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -e error_make_tree_subset_db_proka_r220_qsub.txt
#PBS -o output_make_tree_subset_db_proka_r220_qsub.txt

set -e

GTDB_FOLDER=/rds/general/project/lms-warnecke-raw/live/GTDB_r220
DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220
METADATA="${DB_PROKA_BASE}/db_proka_r220_genus_10s_metadata.csv"

cd $HOME
. load_conda.sh

cd ${HOME}/assembly-main

python -m src.db_proka.make_tree_subset \
	-i $GTDB_FOLDER \
	--metadata_path $METADATA \
    -o $DB_PROKA_BASE \
	> ${HOME}/output_make_tree_subset_db_proka_r220.txt \
	2> ${HOME}/error_make_tree_subset_db_proka_r220.txt
