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

1. Prepare text file with one assembly accession per file

```
GCA_018222585.1
GCF_000337575.1
GCA_016840645.1
GCA_dummy0001.1
```

2. Download the genomes from NCBI

```sh
python -m src.fetch_assemblies -l test_data/assembly_accessions.txt -o test_data
```

Note: output folder specified with argument `-o` must exist.

3. Concatenate all proteins sequences in one fasta files (optional)

```sh
python -m src.postprocessing.concatenate_proteins -i test_data -o test_data/all_proteins.fasta
```

## Outputs

See example output in folder [`test_data/`](test_data), generated with [`run_test.sh`](run_test.sh).
