import os
from pathlib import Path


def get_accession_from_path_name(path : Path):
    return '_'.join(path.name.split('_')[:2])


def get_n_cpus():
    n = os.cpu_count()
    if n is None:
        n = 1
    return n
