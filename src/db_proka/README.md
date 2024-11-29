# DB Proka release 220

DB Proka release 220 (`db_proka_r220`) is a phylogenetically balanced subset of GTDB release 220.

The goal is to build a database that keeps as much of the known diversity as possible while keeping a reasonable balance between the number of genomes in each genus (or another taxonomic level).

## Selection criteria

Selection rules for the main database (`genus_10s`):
- Keep all genera.
- For each genus, select 10 species, one representative genome per species:
  - Type species of genus is always selected,
  - followed by the species with the most genomes (a proxy for "importance").
  - Tie break: CheckM2 completeness (most complete first) and contamination scores (less contaminated first)

The number of species per genus, 10, is a good tradeoff between diversity and phylogenetic balance.

We also provide a smaller database, `family_1s`, that follow the same selection criteria except that the selection happens at family level, and a single species per family is retained.

## Stats

### Main DB (`genus_10s`)

Number of genomes:
- Total    : 68,709 
- Archaea  :  4,649
- Bacteria : 64,060

Number of proteins: 190,445,617

### Small DB (`family_1s`)

Number of genomes:
- Total    : 5,445 
- Archaea  :   564
- Bacteria : 4,896

Number of proteins: 13,888,743

## Files

Files are prefixed with the DB name, i.e. `db_proka_r220_genus_10s_` for the main DB and `db_proka_r220_family_1s_` for the small DB. Main DB files are at the root while the small DB files are in folder `db_proka_r220_family_1s`.

The following files are available in both DB:

- `metadata.csv`
  - Table containing genomes metadata, including taxonomy, genome size, G+C content, etc.
  - Index: `assembly_accession` from NCBI starting with `GCA_` for GenBank and `GCF_` for RefSeq.
- `all_proteins.fasta`
  - All proteins concatenated in a single fasta file.
  - Size (main DB)  : 72 GB
  - Size (small DB) :  5 GB
  - ID format: `<protein_id>@<assembly_accession>$<gtdb_species>`
  - Example: `AAR38856.1@GCA_000008085.1$Nanoarchaeum_equitans`
- `tree_ar53.phyloxml` and `tree_bac120.phyloxml`
  - Subset of archaeal and bacterial representative species tree from GTDB
  - By construction (since we only include representative species) all genomes have an entry.
- `Pfam-A_hits.csv`
  - Concatenated domain hits from Pfam release 37.
  - Columns: `id,assembly_accession,protein_id,hmm_accession,hmm_query,evalue,bitscore,accuracy,start,end`
  - ID has the same format as `all_proteins.fasta`, described above.
  - Size (main DB)  : 39 GB
  - Size (small DB) :  3 GB
- `TIGR_hits.csv`
  - Same as above but for TIGR release 15.
  - Size (main DB)  : 6.8 GB
  - Size (small DB) : 0.5 GB
- `Pfam-A_summary.tsv.gz`
  - Summary table of Pfam domains in all genomes.
  - One row per genome.
  - Index: `assembly_accession` (i.e. first column)
  - Columns: Every single Pfam domains, as per their short name (e.g. `LysM`).
  - Values: Number of proteins where this domain is found in genome.
- `TIGR_summary.tsv.gz`
  - Same as above but for TIGR domains.
- `genomes.tar.gz`
  - All underlying files, including nucleotide sequences and GFF annotations.
  - Size (main DB)  : 398 GB
  - Size (small DB) :  30 GB
  - Genomes folder names within it have the following format: `<assembly_accession>_<assembly_name>`
  - Example: `GCA_000008085.1_ASM808v1` folder content:
    - `GCA_000008085.1_ASM808v1_ani_report.txt`
    - `GCA_000008085.1_ASM808v1_assembly_report.txt`
    - `GCA_000008085.1_ASM808v1_assembly_stats.txt`
    - `GCA_000008085.1_ASM808v1_cds_from_genomic.fna.gz`
    - `GCA_000008085.1_ASM808v1_genomic.fna.gz`
    - `GCA_000008085.1_ASM808v1_genomic.gbff.gz`
    - `GCA_000008085.1_ASM808v1_genomic.gff.gz`
    - `GCA_000008085.1_ASM808v1_Pfam-A.csv.gz`
    - `GCA_000008085.1_ASM808v1_protein.faa.gz`
    - `GCA_000008085.1_ASM808v1_rna_from_genomic.fna.gz`
    - `GCA_000008085.1_ASM808v1_TIGR.csv.gz`


In addition, the main DB contains clustered proteins with 90% sequence identity (computed with `mmseqs linclust`):

- `proteins_clu90.fasta`
  - Cluster representatives
  - 157,980,633 proteins (17% reduction from all proteins)
  - Size: 59 GB
- `proteins_clu90.tsv`
  - Mapping between cluster representatives and members of the cluster.
  - Size: 20 GB
