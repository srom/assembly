"""
Script to download assemlies from [NCBI](https://www.ncbi.nlm.nih.gov/).
"""
import argparse
import logging
import os
import random
import subprocess
import sys
import time
import tempfile
from multiprocessing import Process
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from src.utils import get_n_cpus
from src.ncbi_util.assembly_summary import (
    download_latest_assembly_summary_as_df,
    parse_assembly_summary,
)


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description='Download assemblies from NCBI',
    )
    parser.add_argument(
        '-a', '--assembly', 
        help=(
            'GenBank or RefSeq assembly accession. '
            'See -l to download more than one assembly.'
        ), 
        type=str,
    )
    parser.add_argument(
        '-l', '--assembly_list', 
        help=(
            'Path to file containing one GenbBank or RefSeq '
            'assembly accession per line. '
            'See -a to download a single assembly'
        ), 
        type=str,
    )
    parser.add_argument(
        '-s', '--assembly_summary', 
        help=(
            'Path to NCBI assembly summary file containing'
            'the list of all existing assemblies. '
            'If unspecified, file is downloaded from NCBI\'s ftp.'
        ), 
        type=str,
    )
    parser.add_argument(
        '-o', '--output_folder', 
        help=(
            'Path to folder where output files will be written. '
            'Genomes are be saved to <output_folder>/genomes/<accession>_<asm_name>/'
        ), 
        type=str,
        required=True,
    )
    parser.add_argument(
        '--genbank', 
        help='Convert RefSeq accessions to GenBank accessions.',
        action='store_true', 
        required=False,
    )
    parser.add_argument('--cpu', type=int, default=4)
    args = parser.parse_args()

    assembly = args.assembly
    assemblies_path = args.assembly_list
    summary_path = args.assembly_summary
    output_folder = Path(args.output_folder)
    use_genbank = args.genbank
    n_cpu = min(args.cpu, get_n_cpus())

    if not output_folder.is_dir():
        logger.error(f'Output folder is not a directory: {args.output_folder}')
        sys.exit(1)

    assemblies = get_assembly_list_from_arguments(assembly, assemblies_path)

    if len(assemblies) == 0:
        logger.error('No assemblies to download')
        sys.exit(1)

    if summary_path is None:
        logger.info('No NCBI assembly summary file found - downloading from NCBI (should take a couple of minutes)')
        assembly_summary_df = download_latest_assembly_summary_as_df()
    else:
        logger.info(f'Loading assembly summary from {summary_path}')
        assembly_summary_df = parse_assembly_summary(summary_path)

    download_instructions, missing_accessions = make_download_instructions(assembly_summary_df, assemblies, use_genbank)

    if len(missing_accessions) > 0:
        missing_accessions_path = output_folder / 'missing_accessions.txt'
        with missing_accessions_path.open('w') as f_out:
            for acc in missing_accessions:
                f_out.write(acc + '\n')
        
        logger.warning(f'Assembly accessions not found on NCBI: {len(missing_accessions):,}. Full list stored to {missing_accessions_path}')

    if len(download_instructions) == 0:
        logger.info('No assemblies to be downloaded')
        sys.exit(0)

    logger.info(f'Downloading {len(download_instructions):,} assemblies into {output_folder}')

    n_cpu = min(n_cpu, len(download_instructions))
    n_per_process = int(np.ceil(len(download_instructions) / n_cpu))

    genomes_folder = output_folder / 'genomes'
    genomes_folder.mkdir(exist_ok=True)

    processes = []
    for i in range(n_cpu):
        start = i * n_per_process
        end = start + n_per_process
 
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            download_instructions.iloc[start:end].to_csv(f)
            p = Process(target=worker_main, args=(
                f.name,
                genomes_folder.resolve().as_posix(),
            ))
            p.start()
            processes.append(p)

    for p in processes:
        p.join()

    logger.info('Download completed.')

    metadata_path = output_folder / 'genomes_metadata.csv'
    logger.info(f'Saving genome metadata to {metadata_path}')
    assembly_summary_df.loc[download_instructions.index].to_csv(metadata_path)

    logger.info('DONE')


def worker_main(instructions_path : str, output_folder : str) -> None:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    instructions_path = Path(instructions_path)
    output_folder = Path(output_folder)
    try:
        instructions_df = pd.read_csv(instructions_path, index_col='assembly_accession')

        for i, assembly in enumerate(instructions_df.index):
            instructions = instructions_df.loc[assembly]
            organism_name = instructions.organism_name
            logger.info(f'Downloading assembly {i+1:,} / {len(instructions_df):,}: {assembly} - {organism_name}')
            rsync_url = instructions.ftp_path.replace('https://', 'rsync://')
            try:
                download_assembly_with_retry(rsync_url, output_folder)
            except AssemblyDownloadError as e:
                logger.error(e)
                return
    finally:
        if instructions_path.is_file():
            instructions_path.unlink()


def download_assembly_with_retry(
    rsync_url : str,
    output_folder : os.PathLike,
    n_retries : int = 10,
) -> None:
    try_nb = 0
    while True:
        try_nb += 1
        response = download_asembly(rsync_url, output_folder)

        if response.returncode == 0:
            # Download completed successfully
            break
        
        # Deal with errors:
        # - rsync error code 1 is a user error, raise exception.
        # - otherwise retry with exponential backoff for a set number of times
        stderr_txt = response.stderr.decode('utf-8')
        if response.returncode == 1:
            raise AssemblyDownloadError(stderr_txt) 
        elif try_nb > n_retries:
            raise AssemblyDownloadError(f'{try_nb:,} tries failed, aborting. Error: {stderr_txt}')
        else:
            sleep_seconds = exponential_backoff_sleep_seconds(try_nb)
            logger.warning((
                f'Try number {try_nb:,} failed, retrying in {sleep_seconds:,} seconds. '
                f'Error: {stderr_txt}'
            ))
            time.sleep(sleep_seconds)
            continue


def download_asembly(
    rsync_url : str,
    output_folder : os.PathLike,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            'rsync', 
            '--copy-links', 
            '--recursive',
            '--times', 
            '--quiet',
            '--exclude', 'README.txt',
            rsync_url,
            output_folder.as_posix(),
        ],
        capture_output=True,
    )


def get_assembly_list_from_arguments(
    assembly : Optional[str], 
    assemblies_path : Optional[str]
) -> List[str]:
    if assembly is None and assemblies_path is None:
        logger.error('Specify one of -a or -l to specify assemblies')
        sys.exit(1)
    elif assembly is not None and assemblies_path is not None:
        logger.error('Specify -a or -l but not both')
        sys.exit(1)
    
    if assembly is not None:
        assemblies = [assembly]
    else:
        logger.info(f'Loading assembly list from {assemblies_path}')
        path = Path(assemblies_path)
        if not path.is_file():
            logger.error(f'Fle does not exist: {assemblies_path}')
            sys.exit(1)
        else:
            assemblies = parse_assemblies_from_file(path)

    return assemblies


def parse_assemblies_from_file(path : os.PathLike) -> List[str]:
    with open(path) as f:
        return [
            a.strip() 
            for a in f.readlines()
            if a.strip() != ''
        ]


def make_download_instructions(
    assembly_summary_df : pd.DataFrame, 
    assemblies : List[str],
    use_genbank : bool,
) -> pd.DataFrame:
    """
    Dataframe holding just enough info to be useful for 
    the purpose of downloading genomes.
    """
    genbank_ids = set(assembly_summary_df.index)
    refseq_id_to_genbank = {}
    for genbank_id in assembly_summary_df.index:
        row = assembly_summary_df.loc[genbank_id]
        refseq_id = row['gbrs_paired_asm']
        if pd.isnull(refseq_id) or refseq_id == 'na':
            continue

        refseq_id_to_genbank[refseq_id] = genbank_id

    assembly_ids = set()
    missing_accessions = []
    assembly_use_ref_seq = []
    for assembly in assemblies:
        if assembly in genbank_ids:
            assembly_accession = assembly
        elif assembly in refseq_id_to_genbank:
            assembly_accession = refseq_id_to_genbank[assembly]
            if not use_genbank:
                assembly_use_ref_seq.append((assembly_accession, assembly))
        else:
            missing_accessions.append(assembly)

        assembly_ids.add(assembly_accession)

    # Use RefSeq URL if necessary
    for assembly, refseq_id in assembly_use_ref_seq:
        ftp_path = assembly_summary_df.loc[assembly, 'ftp_path']
        ftp_path = ftp_path.replace(assembly, refseq_id)
        ftp_path = ftp_path.replace('GCA', 'GCF')
        assembly_summary_df.loc[assembly, 'ftp_path'] = ftp_path

    subset_df = assembly_summary_df.loc[sorted(assembly_ids)]

    column_subset = ['organism_name', 'asm_name', 'ftp_path']
    return (
        subset_df[column_subset].copy(),
        missing_accessions,
    )


def exponential_backoff_sleep_seconds(try_nb):
    """
    Return an exponentially larger number of seconds with each try.
    A total of 10 tries adds up to about 5 minutes.
    """
    return np.round(1.6 ** try_nb + random.uniform(0, 1), 1)


class AssemblyDownloadError(Exception):
    pass


if __name__ == '__main__':
    main()
