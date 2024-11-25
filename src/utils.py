import os
from pathlib import Path
import re


def get_accession_from_path_name(path : Path):
    return '_'.join(path.name.split('_')[:2])


def get_n_cpus():
    n = os.cpu_count()
    if n is None:
        n = 1
    return n


def escape_species_name(species_name : str):
    if not isinstance(species_name, str):
        return species_name
    else:
        return re.sub(r'[^a-zA-Z0-9\-_]', '_', species_name.strip())
