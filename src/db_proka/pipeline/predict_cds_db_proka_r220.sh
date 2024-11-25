#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -e error_db_proka_r220_predict_cds_qsub.txt
#PBS -o output_db_proka_r220_predict_cds_qsub.txt

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220

# Input variables can be set with qsub -v
# If not set the following defaults are used:
if [ -z "${BASE_FOLDER}" ]; then
    BASE_FOLDER=$DB_PROKA_BASE
fi
if [ -z "${METADATA}" ]; then
    METADATA="${DB_PROKA_BASE}/db_proka_r220_genus_10s_metadata.csv"
fi

cd $HOME
. load_conda.sh

cd ${HOME}/assembly-main

export PATH="${HOME}/bin:${PATH}"

python -m src.postprocessing.predict_cds \
	-i $BASE_FOLDER \
	--metadata_path $METADATA \
	--cpu 32 \
	> ${HOME}/output_db_proka_r220_predict_cds.txt \
	2> ${HOME}/error_db_proka_r220_predict_cds.txt
