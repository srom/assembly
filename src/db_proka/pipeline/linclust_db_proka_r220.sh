#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=128:mem=128gb
#PBS -e error_linclust_db_proka_r220_qsub.txt
#PBS -o output_linclust_db_proka_r220_qsub.txt

set -e
cd $HOME

DB=/rds/general/project/lms-warnecke-raw/live/db_proka_r220/mmseqs_db/db_proka_r220
DB_clu=/rds/general/project/lms-warnecke-raw/live/db_proka_r220/mmseqs_db/db_proka_r220_clu90
DB_clu_rep=/rds/general/project/lms-warnecke-raw/live/db_proka_r220/mmseqs_db/db_proka_r220_clu90_rep
output_tsv=/rds/general/project/lms-warnecke-raw/live/db_proka_r220/db_proka_r220_genus_10s_proteins_clu90.tsv
output_fasta=/rds/general/project/lms-warnecke-raw/live/db_proka_r220/db_proka_r220_genus_10s_proteins_clu90.fasta
MIN_SEQ_IDENTITY="0.9"

. load_conda.sh
conda activate mmseqs

# Run clustering
mmseqs linclust $DB $DB_clu $TMPDIR \
    --min-seq-id $MIN_SEQ_IDENTITY \
    --threads 128 \
    --split-memory-limit "120G" \
    > output_linclust_db_proka_r220.txt \
    2> error_linclust_db_proka_r220.txt

# Output mapping file
mmseqs createtsv $DB $DB $DB_clu $output_tsv \
    --threads 128 \
    >> output_linclust_db_proka_r220.txt \
    2>> error_linclust_db_proka_r220.txt

# Create representative sub-DB
mmseqs createsubdb $DB_clu $DB $DB_clu_rep \
    >> output_linclust_db_proka_r220.txt \
    2>> error_linclust_db_proka_r220.txt

# Output representative sequences as fasta
mmseqs convert2fasta $DB_clu_rep $output_fasta \
    >> output_linclust_db_proka_r220.txt \
    2>> error_linclust_db_proka_r220.txt
