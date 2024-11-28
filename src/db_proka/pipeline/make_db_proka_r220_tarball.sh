#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -e error_make_db_proka_r220_tarball_qsub.txt
#PBS -o output_make_db_proka_r220_tarball_qsub.txt

set -e
start=$(date +%s)

DB_PROKA_BASE=/rds/general/project/lms-warnecke-raw/live/db_proka_r220

#------ Small DB

SMALL_DB_PATH="${DB_PROKA_BASE}/db_proka_r220_family_1s"
echo "Compressing ${SMALL_DB_PATH}/genomes"
FILENAME="db_proka_r220_family_1s_genomes.tar.gz"
cd $SMALL_DB_PATH

echo "Print the total size of the content in gigabytes before compression"
du -sb genomes/ | awk '{ printf "%.2f GB\n", $1/1024**3 }'

echo "Create tarball of ${SMALL_DB_PATH}/genomes to ${TMPDIR}/${FILENAME}"
tar -czvf - genomes/ > ${TMPDIR}/${FILENAME}

echo "Move ${TMPDIR}/${FILENAME} to ${SMALL_DB_PATH}/${FILENAME}"
mv ${TMPDIR}/${FILENAME} ${SMALL_DB_PATH}/${FILENAME}

chmod 644 ${SMALL_DB_PATH}/${FILENAME}

echo

#------ Full DB

echo "Compressing ${DB_PROKA_BASE}/genomes"
FILENAME="db_proka_r220_genus_10s_genomes.tar.gz"
cd $DB_PROKA_BASE

echo "Print the total size of the content in gigabytes before compression"
du -sb genomes/ | awk '{ printf "%.2f GB\n", $1/1024**3 }'

echo "Create tarball of ${DB_PROKA_BASE}/genomes to ${TMPDIR}/${FILENAME}"
tar -czvf - genomes/ > ${TMPDIR}/${FILENAME}

echo "Move ${TMPDIR}/${FILENAME} to ${DB_PROKA_BASE}/${FILENAME}"
mv ${TMPDIR}/${FILENAME} ${DB_PROKA_BASE}/${FILENAME}

chmod 644 ${DB_PROKA_BASE}/${FILENAME}

#------

end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"
