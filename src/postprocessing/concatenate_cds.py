"""
Concatenate all CDS files in sub folders in one big fasta file.

IDs are changed to "<ID>@<accession>$<species_name>".
(species names are escaped by replacing anything but letters, numbers, hyphen or underscore by an underscore)
"""
import argparse
import logging
import os
import sys
import re
from pathlib import Path
from multiprocessing import Process, Queue
from queue import Empty
import subprocess
import tempfile
from typing import List

import numpy as np
import pandas as pd
from Bio import SeqIO

from src.utils import get_accession_from_path_name, get_n_cpus, escape_species_name


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description='Concatenate CDS fasta files into one big file',
    )
    parser.add_argument(
        '-i', '--base_folder', 
        help='Path to base folder containing genome assemblies, as generated by the script fetch_assemblies.py', 
        type=Path,
        required=True,
    )
    parser.add_argument(
        '-o', '--output_path', 
        help='Path to output file', 
        type=Path,
        required=True,
    )
    parser.add_argument(
        '--metadata_path', 
        help=(
            'Path to GTDB metadata file containing assemblies to be processed. '
            'Defaults to <base_folder>/genomes_metadata.csv'
        ), 
        required=False,
        type=Path,
        default=None,
    )
    parser.add_argument('--cpu', type=int, default=4)
    args = parser.parse_args()

    base_folder = args.base_folder
    output_path = args.output_path
    metadata_path = args.metadata_path
    n_cpus = min(args.cpu, get_n_cpus())

    genomes_folder = base_folder / 'genomes'
    if metadata_path is None:
        metadata_path = base_folder / 'genomes_metadata.csv'

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
            metadata_path,
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
    metadata_path : Path,
):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
        output_path = Path(f.name).resolve()
    
    logger.info(f'Worker {worker_ix+1} is starting')

    metadata_df = pd.read_csv(metadata_path, index_col='assembly_accession')

    for i, path in enumerate(paths):
        if i == 0 or (i+1) % 100 == 0 or (i+1) == len(paths):
            logger.info(f'Worker {worker_ix+1} | Processing assembly {i+1:,} / {len(paths):,}')

        assembly_accession = get_accession_from_path_name(path)
        species_name = metadata_df.loc[assembly_accession, 'gtdb_species']
        species_name_escaped = escape_species_name(species_name)

        cds_path_gz = path / f'{path.name}_cds_from_genomic.fna.gz'
        if not cds_path_gz.is_file():
            logger.error(f'CDS file not found: {cds_path_gz}')
            continue

        # Decompress CDS fasta file
        with tempfile.NamedTemporaryFile(suffix='.fna', delete=False) as cds_file:
            cds_path = Path(cds_file.name).resolve()

            returncode = subprocess.call(
                ['gzip', '-cd', cds_path_gz.resolve().as_posix()],
                stdout=cds_file,
            )
            if returncode != 0:
                logger.error(f'Error while decompressing {cds_path_gz}')
                continue

        try:
            output_records = []
            with cds_path.open() as f_in:
                for record in SeqIO.parse(f_in, 'fasta'):
                    # Add accession and species name ot CDS ID.
                    record.id = f'{record.id}@{assembly_accession}${species_name_escaped}'

                    # Remove any description containing '#' as they may cause issue to downstream processes.
                    # Prodigal CDS descriptions in particular always contain several '#' characters.
                    if record.name is not None and '#' in record.name:
                        record.name = ''
                    if record.description is not None and '#' in record.description:
                        record.description = ''

                    output_records.append(record)

            with output_path.open('a') as f_out:
                SeqIO.write(output_records, f_out, 'fasta')

        finally:
            if cds_path.is_file():
                cds_path.unlink()

    queue.put((worker_ix, output_path))
    logger.info(f'Worker {worker_ix+1} is done')


if __name__ == '__main__':
    main()