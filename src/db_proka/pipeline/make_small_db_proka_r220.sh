#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -e error_make_small_db_proka_r220_qsub.txt
#PBS -o output_make_small_db_proka_r220_qsub.txt

set -e

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220
OUTPUT_FOLDER=/rds/general/project/lms-warnecke-raw/live/db_proka_r220/db_proka_r220_family_1s
METADATA="${OUTPUT_FOLDER}/db_proka_r220_family_1s_metadata.csv"

cd $HOME
. load_conda.sh

cd ${HOME}/assembly-main

python -m src.db_proka.make_small_db \
	-i $DB_PROKA_BASE \
	-m $METADATA \
	-o $OUTPUT_FOLDER \
	> ${HOME}/output_make_small_db_proka_r220.txt \
	2> ${HOME}/error_make_small_db_proka_r220.txt
