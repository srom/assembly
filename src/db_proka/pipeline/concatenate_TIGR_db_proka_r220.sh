#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mem=32gb
#PBS -e error_concatenate_TIGR_db_proka_r220_qsub.txt
#PBS -o output_concatenate_TIGR_db_proka_r220_qsub.txt

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220

# Input variables can be set with qsub -v
# If not set the following defaults are used:
if [ -z "${BASE_FOLDER}" ]; then
    BASE_FOLDER=$DB_PROKA_BASE
fi
if [ -z "${METADATA}" ]; then
    METADATA="${DB_PROKA_BASE}/db_proka_r220_genus_10s_metadata.csv"
fi
if [ -z "${OUTPUT}" ]; then
    OUTPUT="${DB_PROKA_BASE}/db_proka_r220_genus_10s_TIGR_hits.csv"
fi

cd $HOME
. load_conda.sh

cd ${HOME}/assembly-main

python -m src.postprocessing.concatenate_hmm_output \
        -i $BASE_FOLDER \
        -s "TIGR" \
        -o $OUTPUT \
        --metadata_path $METADATA \
        --cpu 32 \
        > ${HOME}/output_concatenate_TIGR_db_proka_r220.txt \
        2> ${HOME}/error_concatenate_TIGR_db_proka_r220.txt
