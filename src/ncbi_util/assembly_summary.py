"""
Utilities to download and parse the assembly summary file from NCBI.

This file contaisn a table of information about all the existing assemblies. 

README: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_assembly_summary.txt
CONTENT: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt 
"""

import os
from pathlib import Path
import tempfile
from typing import Optional

import pandas as pd
import requests

from .urls import assembly_summary_url


def download_latest_assembly_summary_as_df(
    path : Optional[os.PathLike] = None,
    delete_after : bool = True,
) -> pd.DataFrame:
    """
    Download assembly summary from NCBI to the location specified by `path`.
    If `path` is unspecified, download to a temporary file.
    Return parsed file as a pandas DataFrame.
    Optionally removes the downloaded file after DataFrame is created.
    """
    summary_path = None
    try:
        summary_path = download_latest_assembly_summary(path)
        return parse_assembly_summary(summary_path)
    finally:
        if delete_after and summary_path is not None and summary_path.is_file():
            summary_path.unlink()


def download_latest_assembly_summary(
    path : Optional[os.PathLike] = None,
) -> os.PathLike:
    """
    Download assembly summary from NCBI to the location specified by `path`.
    If `path` is unspecified, download to a temporary file.
    Return path to the downloaded file.
    """
    if path is None:
        f = tempfile.NamedTemporaryFile(prefix='ncbi_assembly_summary', delete=False)
        path = Path(f.name).resolve()
    elif isinstance(path, str):
        path = Path(path).resolve()

    with path.open('wb') as f:
        r = requests.get(assembly_summary_url)
        f.write(r.content)

    return path


def parse_assembly_summary(path : os.PathLike) -> pd.DataFrame:
    """
    Parse assembly summary file from NCBI into a pandas DataFrame.
    """
    assembly_summary_df = pd.read_csv(path, sep='\t', skiprows=1)

    columns = assembly_summary_df.columns.tolist()
    columns[0] = 'assembly_accession'
    assembly_summary_df.columns = columns

    return assembly_summary_df.set_index('assembly_accession')
