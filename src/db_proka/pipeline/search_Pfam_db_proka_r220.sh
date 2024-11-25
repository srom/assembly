#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=128:mem=64gb
#PBS -e error_search_Pfam_db_proka_r220_qsub.txt
#PBS -o output_search_Pfam_db_proka_r220_qsub.txt

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220
HMM_DB=/rds/general/project/lms-warnecke-raw/live/Pfam-A_37/Pfam-A.hmm

# Input variables can be set with qsub -v
# If not set the following defaults are used:
if [ -z "${BASE_FOLDER}" ]; then
    BASE_FOLDER=$DB_PROKA_BASE
fi

cd $HOME
. load_conda.sh

cd ${HOME}/assembly-main

python -m src.postprocessing.search_hmm \
	-i $BASE_FOLDER \
	-d $HMM_DB \
	--hmmer_cut_ga \
	--n_processes 120 \
	--n_threads_per_process 1 \
	> ${HOME}/output_search_Pfam_db_proka_r220.txt \
	2> ${HOME}/error_search_Pfam_db_proka_r220.txt
