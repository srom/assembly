"""
Use [HMMER](http://hmmer.org/) or [MMseqs2](https://github.com/soedinglab/mmseqs2) to search protein sequences 
against a database of protein domain profiles (Pfam, TIGRfam, KOfam, etc).
"""
import argparse
import logging
import os
import subprocess
import sys
from multiprocessing import Process
from pathlib import Path
import random
import shutil
import tempfile
from typing import List

import numpy as np
import pandas as pd
from Bio.SearchIO.HmmerIO.hmmer3_domtab import Hmmer3DomtabHmmqueryParser

from src.utils import get_accession_from_path_name


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description=(
            'Use [HMMER](http://hmmer.org/) or [MMseqs2](https://github.com/soedinglab/mmseqs2) to search '
            'protein sequences against a database of protein domain profiles (e.g Pfam, TIGRfam, KOfam, etc). '
            'To search using hmmer, specify either --hmmer_e_value or --cut_ga. '
            'To search using MMseqs2, specify --mmseqs2_sensitivity.'
        ),
    )
    parser.add_argument(
        '-i', '--base_folder', 
        help='Path to base folder containing genome assemblies, as generated by the script fetch_assemblies.py', 
        required=True,
        type=Path,
    )
    parser.add_argument(
        '-d', '--hmm_db', 
        help='Path to profile database', 
        type=Path,
        required=True,
    )
    parser.add_argument(
        '-e', '--hmmer_e_value', 
        type=float,
        default=None,
        help='Domain E-value threshold (hmmer)', 
    )
    parser.add_argument(
        '--hmmer_cut_ga', 
        action='store_true',
        help='Use the GA (gathering) bit scores thresholds instead of E-value (hmmer)', 
    )
    parser.add_argument(
        '--mmseqs2_sensitivity', 
        type=float,
        default=None,
        help='Sensitivity parameter (MMseqs2)', 
    )
    parser.add_argument(
        '-l', '--assembly_subset', 
        help='Path to text file containing one assembly folder name per line', 
        type=str,
        required=False,
    )
    parser.add_argument('--n_processes', type=int, default=2)
    parser.add_argument('--n_threads_per_process', type=int, default=2)
    args = parser.parse_args()

    base_folder = args.base_folder
    hmm_db = args.hmm_db
    hmmer_e_value = args.hmmer_e_value
    hmmer_cut_ga = args.hmmer_cut_ga
    mmseqs2_sensitivity = args.mmseqs2_sensitivity
    n_processes = args.n_processes
    n_threads_per_process = args.n_threads_per_process
    assembly_subset_path = None
    if args.assembly_subset is not None:
        assembly_subset_path = Path(args.assembly_subset)

    if not base_folder.is_dir():
        logger.error(f'Assemblies folder does not exist: {args.assemblies}')
        sys.exit(1)
    
    genomes_folder = base_folder / 'genomes'
    if not genomes_folder.is_dir():
        logger.error(f'Genomes folder does not exist: {genomes_folder}')
        sys.exit(1)
    
    if not hmm_db.is_file():
        logger.error(f'HMM database does not exist: {args.hmm_db}')
        sys.exit(1)

    if not hmmer_cut_ga and hmmer_e_value is None and mmseqs2_sensitivity is None:
        logger.error(f'Specify one of --hmmer_cut_ga, --hmmer_e_value or --mmseqs2_sensitivity')
        sys.exit(1)

    if assembly_subset_path is not None and not assembly_subset_path.is_file():
        logger.error(f'Assembly subset file does not exist: {args.assembly_subset}')
        sys.exit(1)

    assembly_subset = load_assembly_subset(assembly_subset_path)

    use_hmmer = hmmer_cut_ga or hmmer_e_value is not None
    if use_hmmer:
        logger.info(f'Running hmmer: profile db {hmm_db} vs proteome db {base_folder}')
    else:
        logger.info(f'Running MMseqs2: profile db {hmm_db} vs proteome db {base_folder}')

    paths = [
        p for p in genomes_folder.iterdir()
        if (
            p.is_dir() and 
            p.name.startswith('GC') and 
            (
                assembly_subset is None or 
                get_accession_from_path_name(p) in assembly_subset
            )
        )
    ]

    logger.info(f'Total number of assemblies: {len(paths):,}')

    if len(paths) == 0:
        sys.exit(0)

    # Shuffle
    random.shuffle(paths)

    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy HMM DB file to temp directory,
        # i.e. on disk and not on network drive if running on HPC.
        hmm_db_tmp_path = Path(temp_dir) / hmm_db.name
        shutil.copy(hmm_db, hmm_db_tmp_path)

        # Run on tasks on multi processes.
        n_processes = min(n_processes, len(paths))
        n_per_process = int(np.ceil(len(paths) / n_processes))
        processes = []
        for i in range(n_processes):
            start = i * n_per_process
            end = start + n_per_process
    
            p = Process(target=worker_main, args=(
                i,
                hmm_db_tmp_path,
                paths[start:end],
                use_hmmer,
                hmmer_e_value,
                hmmer_cut_ga,
                mmseqs2_sensitivity,
                n_threads_per_process,
            ))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    logger.info('DONE')
    sys.exit(0)


def worker_main(
    worker_ix : int, 
    hmm_db : os.PathLike, 
    paths_to_process : List[os.PathLike], 
    use_hmmer : bool,
    hmmer_e_value : float,
    hmmer_cut_ga : bool,
    mmseqs2_sensitivity : float,
    n_threads_per_process : int,
):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    if use_hmmer:
        hmm_db_name = hmm_db.name.replace('.hmm', '')
    else:
        hmm_db_name = hmm_db.name.replace('.mmseqs2', '')

    for i, path in enumerate(paths_to_process):
        logger.info(f'Worker {worker_ix+1} | Assembly {i+1:,} / {len(paths_to_process)}: {path.name}')

        output_csv_path = path / f'{path.name}_{hmm_db_name}.csv'
        output_csv_path_gz = path / f'{path.name}_{hmm_db_name}.csv.gz'

        if output_csv_path_gz.is_file():
            continue

        if use_hmmer and hmmer_cut_ga:
            threshold_params = ['--cut_ga']
        elif use_hmmer and hmmer_e_value is not None:
            threshold_params = ['--domE', f'{hmmer_e_value}']
        else:
            threshold_params = ['-s', f'{mmseqs2_sensitivity}']

        protein_path_gz = path / f'{path.name}_protein.faa.gz'
        if not protein_path_gz.is_file():
            logger.error(f'Protein file not found: {protein_path_gz}')
            continue

        # Decompress protein fasta file
        with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as protein_file:
            protein_path = Path(protein_file.name).resolve()

            returncode = subprocess.call(
                ['gzip', '-cd', protein_path_gz.resolve().as_posix()],
                stdout=protein_file,
            )
            if returncode != 0:
                logger.error(f'Error while decompressing {protein_path_gz}')
                continue
        
        output_csv_temp_path = Path(tempfile.NamedTemporaryFile(
            suffix=f'_{path.name}_{hmm_db_name}.csv', 
            delete=False,
        ).name)
        domtblout_path = Path(tempfile.NamedTemporaryFile(
            suffix=f'_{path.name}_{hmm_db_name}_domtblout.txt', 
            delete=False,
        ).name)
        try:
            if use_hmmer:
                # Run hmmsearch
                response = run_hmmsearch(
                    hmm_db, 
                    protein_path, 
                    domtblout_path, 
                    threshold_params, 
                    n_threads_per_process,
                )
                if response.returncode != 0:
                    stderr_txt = response.stderr.decode('utf-8')
                    logger.error(f'Error while running `hmmsearch` (hmmer) on {path.name}: {stderr_txt}')
                    continue

                # Process hmmer output into a CSV file
                process_hmmer_output(domtblout_path, output_csv_temp_path)
            else:
                # Run MMseqs2
                response = run_mmseqs2_easy_search(
                    hmm_db, 
                    protein_path, 
                    domtblout_path, 
                    threshold_params, 
                    n_threads_per_process,
                )
                if response.returncode != 0:
                    stderr_txt = response.stderr.decode('utf-8')
                    logger.error(f'Error while running `mmseqs easy-search` on {path.name}: {stderr_txt}')
                    continue

                # Process hmmer output into a CSV file
                process_mmseqs2_output(domtblout_path, output_csv_temp_path)

            # Compress hmmer CSV output to final location
            if output_csv_path.is_file():
                output_csv_path.unlink()
            if output_csv_path_gz.is_file():
                output_csv_path_gz.unlink()

            with output_csv_path_gz.open('wb') as f_out:
                response = subprocess.run(
                    ['gzip', '-c', output_csv_temp_path.resolve().as_posix()], 
                    stdout=f_out,
                )
            if response.returncode != 0:
                logger.error(f'Error while compressing CSV output {output_csv_temp_path}')
                output_csv_path_gz.unlink()
                continue

            # Change permissions
            response = subprocess.run(
                ['chmod', '644', output_csv_path_gz.resolve().as_posix()], 
                capture_output=True,
            )
            if response.returncode != 0:
                stderr_txt = response.stderr.decode('utf-8')
                logger.error(f'Error while setting file permissions for {output_csv_path_gz}: {stderr_txt}')
                continue

        finally:
            if protein_path.is_file():
                protein_path.unlink()
            if output_csv_temp_path.is_file():
                output_csv_temp_path.unlink()
            if domtblout_path.is_file():
                domtblout_path.unlink()


def run_hmmsearch(hmm_db, protein_path, domtblout_path, threshold_params, n_threads_per_process):
    return subprocess.run(
        [
            'hmmsearch',
            '--acc',
            '--noali',
            '-o', '/dev/null',
            '--domtblout', domtblout_path.resolve().as_posix(),
            '--cpu', f'{n_threads_per_process}',
        ] + 
        threshold_params +
        [
            hmm_db.resolve().as_posix(),
            protein_path.resolve().as_posix(),
        ], 
        capture_output=True,
    )


def run_mmseqs2_easy_search(hmm_db, protein_path, domtblout_path, threshold_params, n_threads_per_process):
    tmp_folder = Path(tempfile.gettempdir())
    return subprocess.run(
        [
            'mmseqs',
            'easy-search',
            protein_path.resolve().as_posix(),
            hmm_db.resolve().as_posix(),
            domtblout_path.resolve().as_posix(),
            tmp_folder.resolve().as_posix(),
            '--threads', f'{n_threads_per_process}',
            '--format-output', 'query,target,evalue,bits,qstart,qend',
        ] + 
        threshold_params,
        capture_output=True,
    )


def get_domain_file_columns_dict():
    return {
        'protein_id': [],
        'hmm_accession': [],
        'hmm_query': [],
        'evalue': [],
        'bitscore': [],
        'accuracy': [],
        'start': [],
        'end': [],
    }


def process_hmmer_output(domtblout_path, output_csv_path):
    output_data = get_domain_file_columns_dict()
    with domtblout_path.open() as f:
        parser = Hmmer3DomtabHmmqueryParser(f)

        for record in parser:
            hmm_accession = record.accession
            for protein_id, hit in record.items:
                hmm_query = hit.query_id
                for hit_instance in hit:
                    evalue = hit_instance.evalue
                    bitscore = hit_instance.bitscore
                    start = hit_instance.env_start
                    end = hit_instance.env_end
                    accuracy = hit_instance.acc_avg

                    accession = hmm_accession
                    if accession == '-' or accession is None:
                        accession = hmm_query

                    output_data['protein_id'].append(protein_id)
                    output_data['hmm_accession'].append(accession)
                    output_data['hmm_query'].append(hmm_query)
                    output_data['evalue'].append(evalue)
                    output_data['bitscore'].append(bitscore)
                    output_data['accuracy'].append(accuracy)
                    output_data['start'].append(start)
                    output_data['end'].append(end)

    df = pd.DataFrame.from_dict(output_data).sort_values(['protein_id', 'start'])
    df.to_csv(output_csv_path, index=False)


def process_mmseqs2_output(domtblout_path, output_csv_path):
    output_data = get_domain_file_columns_dict()
    mmseqs2_output_columns = ['query', 'target', 'evalue', 'bitscore', 'start', 'end']
    mmseqs2_output = pd.read_csv(domtblout_path, sep='\t', header=None, names=mmseqs2_output_columns)

    for row in mmseqs2_output.itertuples(index=False):
        output_data['protein_id'].append(row.query)
        output_data['hmm_accession'].append(row.target)
        output_data['hmm_query'].append(row.target)
        output_data['evalue'].append(row.evalue)
        output_data['bitscore'].append(row.bitscore)
        output_data['accuracy'].append(None)
        output_data['start'].append(row.start)
        output_data['end'].append(row.end)

    df = pd.DataFrame.from_dict(output_data).sort_values(['protein_id', 'start'])
    df.to_csv(output_csv_path, index=False)


def load_assembly_subset(assembly_list_path : os.PathLike):
    if assembly_list_path is None:
        return None

    assembly_subset = set()
    with assembly_list_path.open() as f:
        for line in f:
            assembly_folder_name = line.strip()
            if assembly_folder_name != '':
                assembly_subset.add(assembly_folder_name)
    
    return assembly_subset


if __name__ == '__main__':
    main()
