DOC = """
Build a phylogenetically balanced subset of genomes from GTDB r220.
Outputs a list of NCBI accessions ready to be consumed by `fetch_assemblies.py`

GTDB is a fantastic resource, but the depth of some clades makes evolutionary analyses 
prone to strong biases (for instance, there are 39,000 E. coli genomes alone!).

The goal here is to build a database that keeps as much as the known diversity as possible while 
keeping a reasonable balance between the number of genomes in each genus.

Selection rules:
- Keep all genera.
- For each genus, select N=10 species, one representative genome per species:
  - Type species of genus is always selected,
  - followed by the species with the most genomes (a proxy for "importance").
  - Tie break: CheckM2 completeness (most complete first) and contamination scores (less contaminated first)

The number of species per genus (N) is configurable but 10 is a good tradeoff between diversity and phylogenetic balance.

NCBI accession versions may have changed since the release of GTDB, therefore this script checks 
if a newer version is available from the NCBI master list and removes outdated accessions.
"""
import argparse
import logging
from pathlib import Path
import requests
import sys
import tempfile

import pandas as pd

from src.ncbi_util.assembly_summary import parse_assembly_summary, download_latest_assembly_summary_as_df


# High level constants
GTDB_VERSION = '220.0'
GTDB_BASE_PATH = 'https://data.ace.uq.edu.au/public/gtdb/data/releases'
DEFAULT_N_SPECIES_PER_GENUS = 10

# Transformed constants
GTDB_VERSION_SHORT = GTDB_VERSION.split('.')[0]
GTDB_METADATA_ARCHAEA = GTDB_BASE_PATH + (
    f'/release{GTDB_VERSION_SHORT}/{GTDB_VERSION}/ar53_metadata_r{GTDB_VERSION_SHORT}.tsv.gz'
)
GTDB_METADATA_BACTERIA = GTDB_BASE_PATH + (
    f'/release{GTDB_VERSION_SHORT}/{GTDB_VERSION}/bac120_metadata_r{GTDB_VERSION_SHORT}.tsv.gz'
)

logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(description=DOC, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-o', '--output_folder', 
        type=Path,
        required=True,
    )
    parser.add_argument(
        '--gtdb_metadata',
        help=f'Path to GTDB metadata file. Version {GTDB_VERSION} gets downloaded if no file is provided.',
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        '-s', '--assembly_summary', 
        help=(
            'Path to NCBI assembly summary file containing'
            'the list of all existing assemblies. '
            'If unspecified, file is downloaded from NCBI\'s ftp.'
        ), 
        type=Path,
        required=False,
        default=None,
    )
    parser.add_argument(
        '-n', '--n_species_per_genus', 
        type=int,
        required=False,
        default=DEFAULT_N_SPECIES_PER_GENUS,
    )
    args = parser.parse_args()

    output_folder = args.output_folder
    gtdb_metadata_path = args.gtdb_metadata
    assembly_summary_path = args.assembly_summary
    n_species_per_genus = args.n_species_per_genus

    if not output_folder.is_dir():
        logger.error(f'Output folder does not exist: {args.output_folder}')
        sys.exit(1)
    elif gtdb_metadata_path is not None and not gtdb_metadata_path.is_file():
        logger.error(f'GTDB metadata file does not exist: {args.gtdb_metadata}')
        sys.exit(1)
    elif assembly_summary_path is not None and not assembly_summary_path.is_file():
        logger.error(f'NCBI assembly summary file does not exist: {args.assembly_summary}')
        sys.exit(1)

    logger.info('Building a phylogenetically balanced subset of genomes from GTDB r220')
    logger.info(f'output_folder       = {output_folder}')
    logger.info(f'gtdb_metadata_path  = {gtdb_metadata_path}')
    logger.info(f'n_species_per_genus = {n_species_per_genus:,}')

    gtdb_metadata = load_or_download_gtdb_metadata(gtdb_metadata_path, output_folder)

    if assembly_summary_path is None:
        logger.info('No NCBI assembly summary file found - downloading from NCBI (should take a couple of minutes)')
        assembly_summary_df = download_latest_assembly_summary_as_df(
            output_folder / 'ncbi_assembly_summary.txt',
            delete_after=False,
        )
    else:
        logger.info(f'Loading assembly summary from {assembly_summary_path}')
        assembly_summary_df = parse_assembly_summary(assembly_summary_path)

    logger.info('Selecting genomes')
    accessions = make_selection(gtdb_metadata, n_species_per_genus)
    gtdb_subset = gtdb_metadata.loc[accessions]

    logger.info('Updating accession versions')
    gtdb_subset = update_accessions_with_ncbi_metadata(gtdb_subset, assembly_summary_df)

    # Check that we indeed only have a single genome per species.
    assert len(gtdb_subset.drop_duplicates('gtdb_species')) == len(gtdb_subset)

    logger.info('Outputting list of accessions and metadata')
    gtdb_subset.reset_index()[['assembly_accession']].to_csv(
        output_folder / f'gtdb_metadata_r{GTDB_VERSION_SHORT}_subset_accessions.txt', 
        index=False,
        header=False,
    )
    gtdb_subset.to_csv(output_folder / f'gtdb_metadata_r{GTDB_VERSION_SHORT}_subset.csv')

    logger.info('Stats:')
    n_archaea = len(gtdb_subset[gtdb_subset['domain'] == 'Archaea'])
    n_bacteria = len(gtdb_subset[gtdb_subset['domain'] == 'Bacteria'])
    total = n_archaea + n_bacteria
    assert len(gtdb_subset) == total
    logger.info(f'\t# Archaea  : {n_archaea:,}')
    logger.info(f'\t# Bacteria : {n_bacteria:,}')
    logger.info(f'\t# Genomes  : {total:,}')

    logger.info('DONE')
    sys.exit(0)


def make_selection(gtdb_metadata_input, n_species_per_genus):
    accession_list = []

    # Add number of genomes per species
    gtdb_metadata = gtdb_metadata_input[[
        'gtdb_genus', 
        'gtdb_species',
        'gtdb_type_species_of_genus', 
        'gtdb_representative',
        'checkm2_completeness',
        'checkm2_contamination',
    ]].copy()
    species_stats = gtdb_metadata.reset_index()[
        ['gtdb_species', 'assembly_accession']
    ].groupby('gtdb_species').count()
    n_genomes_in_species = {
        species: species_stats.loc[species]['assembly_accession'] for species in species_stats.index
    }
    gtdb_metadata['n_genomes_in_species'] = gtdb_metadata['gtdb_species'].apply(lambda s: n_genomes_in_species[s])

    # Sort metadata appropriately
    gtdb_metadata_sorted = gtdb_metadata[
        # Species representative only
        gtdb_metadata['gtdb_representative'] == 't'
    ].sort_values(
        [
            'gtdb_genus', 
            'gtdb_type_species_of_genus', 
            'n_genomes_in_species',
            'checkm2_completeness',
            'checkm2_contamination',
        ], 
        ascending=[
            True,   # Sort by genus
            False,  # Genus type strain first
            False,  # Higher number of genomes first
            False,  # tiebreak 1: most complete genomes first
            True,   # tiebreak 2: least amount of contamination first
        ]
    ).reset_index().set_index('gtdb_genus', drop=True)

    # Select top N genomes per genus
    genera = sorted(set(gtdb_metadata_sorted.index))
    for genus in genera:
        selection_df = gtdb_metadata_sorted.loc[[genus]].head(n_species_per_genus)

        accession_list.extend(
            sorted(selection_df['assembly_accession'].unique())
        )

    return accession_list


def load_or_download_gtdb_metadata(gtdb_metadata_path : Path, output_folder : Path) -> pd.DataFrame:
    if gtdb_metadata_path is not None:
        logger.info(f'Loading metadata from {gtdb_metadata_path} ')
        return pd.read_csv(gtdb_metadata_path, index_col='assembly_accession')
    else:
        return download_gtdb_metadata(output_folder)
    

def download_gtdb_metadata(output_folder : Path) -> pd.DataFrame:
    output_path = output_folder / f'gtdb_r{GTDB_VERSION_SHORT}_metadata_all.csv.gz'
    logger.info(f'Downloading GTDB metadata version {GTDB_VERSION} to {output_path}')

    archaea_df = download_tsv_from_gtdb(is_archaea=True)
    bacteria_df = download_tsv_from_gtdb(is_archaea=False)

    metadata_df = pd.concat([archaea_df, bacteria_df], ignore_index=True)

    metadata_df = metadata_df.rename(columns={'accession': 'gtdb_accession'})
    metadata_df['assembly_accession'] = metadata_df['gtdb_accession'].apply(lambda a: a[3:])
    metadata_df = metadata_df.set_index('assembly_accession', drop=True)

    metadata_df['domain'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 'd__'))
    metadata_df['gtdb_phylum'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 'p__'))
    metadata_df['gtdb_class'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 'c__'))
    metadata_df['gtdb_order'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 'o__'))
    metadata_df['gtdb_family'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 'f__'))
    metadata_df['gtdb_genus'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 'g__'))
    metadata_df['gtdb_species'] = metadata_df['gtdb_taxonomy'].apply(lambda t: extract_taxonomy(t, 's__'))

    metadata_df.to_csv(output_path, compression='gzip')

    return metadata_df


def download_tsv_from_gtdb(is_archaea : bool):
    url = GTDB_METADATA_ARCHAEA if is_archaea else GTDB_METADATA_BACTERIA

    f = tempfile.NamedTemporaryFile(prefix='gtdb_metadata', delete=False)
    path = Path(f.name).resolve()

    try:
        with path.open('wb') as f:
            r = requests.get(url)
            f.write(r.content)
        
        return pd.read_csv(path, sep='\t', compression='gzip')
    
    finally:
        if path.is_file():
            path.unlink()


def extract_taxonomy(tax : str, level : str):
    for t in tax.split(';'):
        if t.startswith(level):
            return t.replace(level, '')


def update_accessions_with_ncbi_metadata(gtdb_subset, assembly_summary_df):
    """
    Since GTDB processed the data, the underlying NCBI genomes may have changed in one of two ways:
    - update from GCA_XXXXX.1 to GCA_XXXXX.2 (version bump)
    - deletion (genome no longer available on NCBI)

    This function deals with both possibilities and update the metdata subset file accordingly.
    """
    accessions_from_ncbi = set(assembly_summary_df.index)
    accessions_from_ncbi_no_version = {
        a.split('.')[0]
        for a in accessions_from_ncbi
    }

    accessions_not_in_ncbi = sorted(set(gtdb_subset.index) - accessions_from_ncbi)
    accessions_not_in_ncbi_no_version = {
        a.split('.')[0]
        for a in accessions_not_in_ncbi
    }
    version_change_no_version = accessions_not_in_ncbi_no_version & accessions_from_ncbi_no_version

    assembly_summary_df['accession_no_version'] = [a.split('.')[0] for a in assembly_summary_df.index]
    assembly_summary_reindexed = assembly_summary_df.reset_index().set_index('accession_no_version', drop=True)

    not_available = {
        a for a in sorted(accessions_not_in_ncbi_no_version) 
        if (
            a in accessions_not_in_ncbi and
            a.split('.')[0] not in version_change_no_version
        )
    }
    version_change = {
        old_acc: assembly_summary_reindexed.loc[old_acc.split('.')[0], 'assembly_accession']
        for old_acc in gtdb_subset.index
        if (
            old_acc in accessions_not_in_ncbi and 
            old_acc.split('.')[0] in version_change_no_version
        )
    }

    logger.info(f'Number of accessions with a version bump: {len(version_change):,}')
    logger.info(f'Number of deleted accessions: {len(not_available):,}')

    if len(version_change) > 0 or len(not_available) > 0:
        accessions_to_keep = [a for a in gtdb_subset.index if a not in not_available]
        gtdb_subset_copy = gtdb_subset.loc[accessions_to_keep]

        gtdb_subset_copy.index = pd.Index(
            [
                version_change[a] if a in version_change else a
                for a in gtdb_subset_copy.index
            ],
            name=gtdb_subset_copy.index.name,
        )
        return gtdb_subset_copy
    else:
        return gtdb_subset


if __name__ == '__main__':
    main()
