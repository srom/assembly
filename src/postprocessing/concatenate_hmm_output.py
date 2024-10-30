"""
Concatenate CSV files recapitulating Pfam or TIGR hits in each genome.
"""
import argparse
import logging
import os
import sys
from pathlib import Path
from multiprocessing import Process, Queue
from queue import Empty
import subprocess
import tempfile
from typing import List

import numpy as np
import pandas as pd

from src.utils import get_accession_from_path_name, get_n_cpus


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description='Concatenate CSV files recapitulating Pfam or TIGR hits in each genome.',
    )
    parser.add_argument(
        '-i', '--assemblies', 
        help='Path to base folder containing assemblies', 
        type=str,
        required=True,
    )
    parser.add_argument(
        '-s', '--suffix',
        help='Name of the source', 
        choices=['Pfam-A', 'TIGR'],
        type=str,
        required=True,
    )
    parser.add_argument(
        '-o', '--output_path', 
        help='Path to output file', 
        type=str,
        required=True,
    )
    parser.add_argument('--cpu', type=int, default=4)
    args = parser.parse_args()

    base_folder = Path(args.assemblies)
    suffix = Path(args.suffix)
    output_path = Path(args.output_path)
    n_cpus = min(args.cpu, get_n_cpus())

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
    elif output_path.is_file():
        logger.error(f'Output file already exists: {args.output_path}')
        sys.exit(1)

    logger.info('Loading metadata')
    metadata_df = pd.read_csv(metadata_path, index_col='assembly_accession')
    accessions = set(metadata_df.index)

    logger.info('Finding relevant assembly folders')
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

    logger.info(f'Concatenating {suffix} files for each genomes into one file at {output_path}')

    n_processes = min(n_cpus, len(paths))
    n_per_process = int(np.ceil(len(paths) / n_processes))

    processes = []
    queue = Queue()
    for i in range(n_processes):
        start = i * n_per_process
        end = start + n_per_process
 
        p = Process(target=worker_main, args=(
            i,
            suffix,
            paths[start:end],
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

    # Keep original ordering (by accession)
    partial_paths = sorted(partial_paths, key=lambda t: t[0])

    # Concatenate files
    try:
        with output_path.open('w') as f_out:
            sorted_paths = [p.as_posix() for _, p in partial_paths]
            returncode = subprocess.call(
                ['cat'] + sorted_paths, 
                stdout=f_out,
            )
            if returncode != 0:
                logger.error(f'Error while concatenating files')
                sys.exit(1)

    finally:
        for _, p in partial_paths:
            if p.is_file():
                p.unlink()

    # Change permissions
    response = subprocess.run(
        ['chmod', '644', output_path.resolve().as_posix()], 
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
    suffix : str,
    paths : List[os.PathLike], 
    queue : Queue,
):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
        output_path = Path(f.name).resolve()
    
    logger.info(f'Worker {worker_ix+1} is starting')

    include_header = worker_ix == 0

    for i, path in enumerate(paths):
        if i == 0 or (i+1) % 100 == 0 or (i+1) == len(paths):
            logger.info(f'Worker {worker_ix+1} | Processing assembly {i+1:,} / {len(paths):,}')

        assembly_accession = get_accession_from_path_name(path)

        hmm_output_path_gz = path / f'{path.name}_{suffix}.csv.gz'
        if not hmm_output_path_gz.is_file():
            logger.error(f'File not found: {hmm_output_path_gz}')
            continue

        # Decompress protein fasta file
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as hmm_output_file:
            hmm_output_path = Path(hmm_output_file.name).resolve()

            returncode = subprocess.call(
                ['gzip', '-cd', hmm_output_path_gz.resolve().as_posix()],
                stdout=hmm_output_file,
            )
            if returncode != 0:
                logger.error(f'Error while decompressing {hmm_output_path_gz}')
                continue

        try:
            df = pd.read_csv(hmm_output_path)
            columns = df.columns.tolist()
            df['id'] = df['protein_id'].apply(lambda protein_id: f'{protein_id}@{assembly_accession}')
            df['assembly_accession'] = assembly_accession
            columns = ['id', 'assembly_accession'] + columns
            df[columns].to_csv(output_path, index=False, header=include_header)

        finally:
            if hmm_output_path.is_file():
                hmm_output_path.unlink()

    queue.put((worker_ix, output_path))
    logger.info(f'Worker {worker_ix+1} is done')


if __name__ == '__main__':
    main()
