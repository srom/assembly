Assembly
----------

Fetch and process genome assembly data from [NCBI](https://www.ncbi.nlm.nih.gov/).


## Installation

### Locally with conda

```sh
git clone https://github.com/srom/assembly.git
cd assembly
conda env create -f environment.yml
conda activate assembly
```

### Pip

PyPI package TBD.

### Conda

Conda package TBD.

## Run

### Download genomes from NCBI

Prepare text file with one assembly accession per file, e.g.:

```
GCA_018222585.1
GCF_000337575.1
GCA_016840645.1
GCA_dummy0001.1
```

Non-existent accessions (such as the last one above) are logged and reported in a file (`missing-accessions.txt`), but they won't make the execution fail and other genomes will still get downloaded. 

Run command:

```sh
python -m src.fetch_assemblies -l test_data/assembly_accessions.txt -o test_data
```

### Predict coding sequences (CDS) with Prodigal (optional)

Automatically only runs on genomes without protein fasta file available. Suitable for prokaryotes or phages. [Prodigal](https://github.com/hyattpd/Prodigal) must be installed.

```sh
python -m src.postprocessing.predict_cds -i test_data
```

### Concatenate all proteins sequences in one fasta file (optional)

```sh
python -m src.postprocessing.concatenate_proteins -i test_data -o test_data/all_proteins.fasta
```

## Outputs

See example output in folder [`test_data/`](test_data), generated with [`run_test.sh`](run_test.sh).
