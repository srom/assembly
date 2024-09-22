"""
Concatenate all protein files in sub folders in one big fasta file.

Protein IDs are changed from "<protein_id>" to "<protein_id>@<accession>".
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
from Bio import SeqIO

from src.utils import get_accession_from_path_name, get_n_cpus


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description='Concatenate protein fasta files into one big file',
    )
    parser.add_argument(
        '-i', '--assemblies', 
        help='Path to base folder containing assemblies', 
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
    output_path = Path(args.output_path)
    n_cpus = min(args.cpu, get_n_cpus())

    metadata_path = base_folder / 'genomes_metadata.csv'

    if not base_folder.is_dir():
        logger.error(f'Assemblies folder does not exist: {args.assemblies}')
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
        p for p in base_folder.iterdir()
        if (
            p.is_dir() and 
            p.name.startswith('GC') and 
            get_accession_from_path_name(p) in accessions
        )
    ])

    logger.info(f'Total number of assemblies: {len(paths):,}')

    if len(paths) == 0:
        sys.exit(0)

    n_processes = min(n_cpus, len(paths))
    n_per_process = int(np.ceil(len(paths) / n_processes))

    processes = []
    queue = Queue()
    for i in range(n_processes):
        start = i * n_per_process
        end = start + n_per_process
 
        p = Process(target=worker_main, args=(
            i,
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
    paths : List[os.PathLike], 
    queue : Queue,
):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
        output_path = Path(f.name).resolve()
    
    logger.info(f'Worker {worker_ix+1} is starting')

    for i, path in enumerate(paths):
        if i == 0 or (i+1) % 100 == 0 or (i+1) == len(paths):
            logger.info(f'Worker {worker_ix+1} | Processing assembly {i+1:,} / {len(paths):,}')

        assembly_accession = get_accession_from_path_name(path)

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

        try:
            output_records = []
            with protein_path.open() as f_in:
                for record in SeqIO.parse(f_in, 'fasta'):
                    record.id = f'{record.id}@{assembly_accession}'
                    output_records.append(record)

            with output_path.open('a') as f_out:
                SeqIO.write(output_records, f_out, 'fasta')

        finally:
            if protein_path.is_file():
                protein_path.unlink()

    queue.put((worker_ix, output_path))
    logger.info(f'Worker {worker_ix+1} is done')


if __name__ == '__main__':
    main()
