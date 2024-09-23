"""
Summarize domain hits created with script `search_hmm.py` into a single file 
with one row per assembly and one column per domain containing the hit count.

Columns:
assembly_accession = GCA_XXX
asm_name = project_xyz
domain_1 = 0
domain_2 = 3
domain_3 = 1
...
domain_N = 0 
"""
import argparse
import logging
import os
import re
import sys
from pathlib import Path
from multiprocessing import Process, Queue
from queue import Empty
import subprocess
import tempfile
from typing import Dict, List, Union

import numpy as np
import pandas as pd

from src.utils import get_accession_from_path_name


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description='Count protein domains in assemblies',
    )
    parser.add_argument(
        '-i', '--assemblies', 
        help='Path to base folder containing assemblies', 
        type=str,
        required=True,
    )
    parser.add_argument(
        '-d', '--hmm_db', 
        help='Path to HMM database', 
        type=str,
        required=True,
    )
    parser.add_argument(
        '-o', '--output_path', 
        help='Path to output file', 
        type=str,
        required=True,
    )
    parser.add_argument('--n_cpus', type=int, default=4)
    args = parser.parse_args()

    base_folder = Path(args.assemblies)
    hmm_db = Path(args.hmm_db)
    output_path = Path(args.output_path)
    n_cpus = args.n_cpus

    metadata_path = base_folder / 'genomes_metadata.csv'
    genomes_folder = base_folder / 'genomes'

    if not base_folder.is_dir():
        logger.error(f'Assemblies folder does not exist: {args.assemblies}')
        sys.exit(1)
    elif not genomes_folder.is_dir():
        logger.error(f'Genomes folder does not exist: {genomes_folder}')
        sys.exit(1)
    elif not metadata_path.is_file():
        logger.error(f'Metadata file does not exist: {metadata_path}')
        sys.exit(1)
    elif not hmm_db.is_file():
        logger.error(f'HMM database does not exist: {args.hmm_db}')
        sys.exit(1)
    elif output_path.is_file():
        logger.error(f'Output file already exists: {args.output_path}')
        sys.exit(1)

    logger.info('Loading metadata')
    metadata_df = pd.read_csv(metadata_path, index_col='assembly_accession')
    accessions = set(metadata_df.index)

    hmm_db_name = hmm_db.name.replace('.hmm', '')
    logger.info(f'Counting domains: HMM db {hmm_db_name} vs protein db @ {base_folder}')

    paths = sorted([
        p for p in genomes_folder.iterdir()
        if (
            p.is_dir() and 
            p.name.startswith('GC') and 
            get_accession_from_path_name(p) in accessions
        )
    ])

    logger.info(f'Total number of assemblies: {len(paths):,}')
    if len(paths) == 0:
        sys.exit(0)

    domain_columns = get_domain_columns(hmm_db)

    n_processes = min(n_cpus, len(paths))
    n_per_process = int(np.ceil(len(paths) / n_processes))

    processes = []
    queue = Queue()
    for i in range(n_processes):
        start = i * n_per_process
        end = start + n_per_process
 
        p = Process(target=worker_main, args=(
            i,
            hmm_db,
            paths[start:end],
            domain_columns,
            queue,
        ))
        p.start()
        processes.append(p)

    partial_paths = []
    for p in processes:
        p.join()
        try:
            temp_path = queue.get_nowait()
            partial_paths.append(temp_path)
        except Empty:
            continue

    # Sort path by worker index (because only first worker contains the csv header)
    partial_paths = sorted(partial_paths, key=lambda t: t[0])

    # Concatenate files
    try:
        with output_path.open('w') as f:
            sorted_paths = [p.as_posix() for _, p in partial_paths]
            returncode = subprocess.call(
                ['cat'] + sorted_paths, 
                stdout=f,
            )
            if returncode != 0:
                logger.error(f'Error while concatenating files')
                sys.exit(1)

    finally:
        for _, p in partial_paths:
            if p.is_file():
                p.unlink()

    # Compress output file
    response = subprocess.run(['gzip', output_path.resolve().as_posix()], capture_output=True)
    if response.returncode != 0:
        stderr_txt = response.stderr.decode('utf-8')
        logger.error(f'Error while compressing CSV output {output_path}: {stderr_txt}')
        sys.exit(1)

    # Change permissions
    response = subprocess.run(
        ['chmod', '644', f'{output_path.resolve().as_posix()}.gz'], 
        capture_output=True,
    )
    if response.returncode != 0:
        stderr_txt = response.stderr.decode('utf-8')
        logger.error(f'Error while setting file permissions: {stderr_txt}')
        sys.exit(1)

    logger.info('DONE')
    sys.exit(0)


def worker_main(
    worker_ix : int, 
    hmm_db : os.PathLike, 
    paths : List[os.PathLike], 
    domain_columns : List[str], 
    queue : Queue,
):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as f:
        output_path = Path(f.name).resolve()

    header = worker_ix == 0
    
    logger.info(f'Worker {worker_ix+1} starting')
    hmm_db_name = hmm_db.name.replace('.hmm', '')
    data = init_data_holder(domain_columns)
    save_every = 100
    for i, path in enumerate(paths):
        if i == 0 or (i+1) % save_every == 0 or (i+1) == len(paths):
            logger.info(f'Worker {worker_ix + 1} | Processing assembly {i+1:,} / {len(paths):,}')

        assembly_accession, asm_name = extract_metadata_from_path_name(path.name)

        data['assembly_accession'].append(assembly_accession)
        data['asm_name'].append(asm_name)

        domains_path = path / f'{path.name}_{hmm_db_name}.csv.gz'
        domains_df = pd.read_csv(
            domains_path, 
            usecols=['hmm_query', 'protein_id'], 
        )

        # Keep one domain hit per protein
        domains_df = domains_df.drop_duplicates(['protein_id', 'hmm_query']).set_index('hmm_query')

        domains_ix = domains_df.index
        for domain in domain_columns:
            count = (domains_ix == domain).sum()
            data[domain].append(count)

        n_records = len(data['assembly_accession'])
        if n_records >= save_every:
            append_to_output_file(data, output_path, header)
            data = init_data_holder(domain_columns)
            header = False

    n_records = len(data['assembly_accession'])
    if n_records > 0:
        append_to_output_file(data, output_path, header)

    queue.put((worker_ix, output_path))


def append_to_output_file(data : Dict[str, List[Union[str, int]]], output_path : os.PathLike, header : bool):
    df = pd.DataFrame.from_dict(data)
    df.to_csv(
        output_path, 
        sep='\t', 
        mode='a', 
        header=header, 
        index=False,
    )


def init_data_holder(domain_columns : List[str]):
    data = {
        'assembly_accession': [],
        'asm_name': [],
    }
    for domain_column in domain_columns:
        data[domain_column] = []

    return data


def get_domain_columns(hmm_db : os.PathLike):
    with tempfile.NamedTemporaryFile(suffix='.hmmstat.txt', delete=False) as f:
        stats_path = Path(f.name).resolve()

        returncode = subprocess.call(
            [
                'hmmstat',
                hmm_db.resolve().as_posix(),
            ],
            stdout=f,
        )
        if returncode != 0:
            raise ValueError(f'Error while running hmmstat')

    try:
        df = pd.read_csv(
            stats_path, 
            sep='\s+', 
            comment='#',
            header=None,
            names=[
                'idx', 'name', 'accession', 'nseq', 'eff_nseq', 
                'M', 'relent', 'info', 'p relE', 'compKL',
            ],
            usecols=['name', 'accession'],
            index_col='name',
        )
        return sorted(df.index.tolist())

    finally:
        if stats_path.is_file():
            stats_path.unlink()


def extract_metadata_from_path_name(path_name : str):
    m = re.match(r'^(GC[AF]_[^_]+)_(.+)$', path_name)
    return m[1], m[2]


if __name__ == '__main__':
    main()
