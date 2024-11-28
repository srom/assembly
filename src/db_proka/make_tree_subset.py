"""
Create a subset of the GTDB trees of species representative (one for archaea, one for bacteria).

The subset only includes genomes included in our DB. The prefix RS_ or GB_ is scrapped from leaf IDs.
"""
import argparse
import copy
import logging
from pathlib import Path
import sys
from typing import Set

import pandas as pd
from Bio import Phylo


logger = logging.getLogger()


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(processName)-10s (%(levelname)s) %(message)s')

    parser = argparse.ArgumentParser(
        description='Create subset of GTDB trees',
    )
    parser.add_argument(
        '-i', '--gtdb_folder', 
        help='Path to GTDB folder containing raw trees', 
        type=Path,
        required=True,
    )
    parser.add_argument(
        '-m', '--metadata_path', 
        help='Path to GTDB metadata file containing assemblies to be processed.', 
        required=True,
        type=Path,
    )
    parser.add_argument(
        '-o', '--output_folder', 
        help='Path to output folder', 
        type=Path,
        required=True,
    )
    args = parser.parse_args()

    gtdb_folder = args.gtdb_folder
    metadata_path = args.metadata_path
    output_folder = args.output_folder

    if not gtdb_folder.is_dir():
        logger.error(f'GTDB folder does not exist: {args.gtdb_folder}')
        sys.exit(1)
    elif not metadata_path.is_file():
        logger.error(f'Metadata file does not exist: {metadata_path}')
        sys.exit(1)
    if not output_folder.is_dir():
        logger.error(f'Output folder does not exist: {args.output_folder}')
        sys.exit(1)
    
    archaeal_tree_path = gtdb_folder / 'ar53_r220.tree'
    bacterial_tree_path = gtdb_folder / 'bac120_r220.tree'

    if not archaeal_tree_path.is_file():
        logger.error(f'Archaeal tree does not exist: {archaeal_tree_path}')
        sys.exit(1)
    elif not bacterial_tree_path.is_file():
        logger.error(f'Bacterial tree does not exist: {bacterial_tree_path}')
        sys.exit(1)

    logger.info('Loading metadata')
    metadata_df = pd.read_csv(metadata_path, index_col='assembly_accession')
    gtdb_accessions = set(metadata_df['gtdb_accession'].unique())

    logger.info(f'Loading trees from {gtdb_folder}')
    arc_tree = Phylo.read(archaeal_tree_path, 'newick')
    bac_tree = Phylo.read(bacterial_tree_path, 'newick')

    for (domain, tree) in [('archaea', arc_tree), ('bacteria', bac_tree)]:
        output_path = output_folder / f'tree_{domain}.phyloxml'
        logger.info(f'Processing {domain} tree')

        updated_tree = prune_leaves_with_unknown_id(tree, gtdb_accessions)
        remove_gtdb_prefix_from_leaf_ids(updated_tree)

        n_genomes = len(list(updated_tree.get_terminals()))
        logger.info(f'Number of genomes ({domain}): {n_genomes:,}')

        logger.info(f'Exporting to {output_path}')
        with output_path.open('w') as f_out:
            Phylo.write([updated_tree], f_out, 'phyloxml')

    logger.info('DONE')
    sys.exit(0)


def prune_leaves_with_unknown_id(
    tree: Phylo.BaseTree.Tree, 
    id_set: set[str],
) -> Phylo.BaseTree.Tree:
    """
    Prune tree leaves that do not match any of the known ids.
    Returns a copy of the original tree and collapses single-child clades.
    """
    def prune_and_collapse(clade, parent_clade=None):
        """
        Recursively prune clades that don't match the id_set and collapse single-child clades.
        """
        # Iterate over a copy of the clades list to avoid modifying it while iterating
        for child in clade.clades[:]:
            prune_and_collapse(child, clade)
        
        # If this clade has no children after pruning, and it's not in id_set, prune it
        if not clade.clades and clade.name not in id_set:
            if parent_clade:
                parent_clade.clades.remove(clade)

        # If the clade has only one child, collapse it into that child
        if len(clade.clades) == 1:
            only_child = clade.clades[0]
            clade.branch_length = (clade.branch_length or 0) + (only_child.branch_length or 0)
            clade.clades = only_child.clades
            clade.name = only_child.name

    # Deep copy the tree to avoid modifying the original
    out_tree = copy.deepcopy(tree)

    # Start pruning and collapsing from the root
    prune_and_collapse(out_tree.root)

    return out_tree


def remove_gtdb_prefix_from_leaf_ids(tree):
    for leaf in tree.get_terminals():
        leaf.name = leaf.name[3:]


if __name__ == '__main__':
    main()
