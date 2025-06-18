
# Primer Design Pipeline


This project implements an automated pipeline for primer design from public genomic sequences, designed for DNA barcoding studies. The workflow is built with [**Snakemake**](https://snakemake.readthedocs.io/en/stable/) and integrates tools such as [**Biopython**](https://biopython.org/wiki/Documentation), Clustal Omega, EMBOSS, and Primer3.

## Requeriments and specifications

This pipeline depends on a `virtualenv` Python environment and a GNU/Linux system (specifically developed on Ubuntu 24.04.1).  
External command-line tools are required and must be accessible via the system `PATH`:

| Package                                                       | Command called in Snakefile | Version |
| ------------------------------------------------------------- | --------------------------- | ------- |
| [Clustal Omega ](http://www.clustal.org/omega/clustalo-api/)  | `clustalo`                  | 1.2.4   |
| [EMBOSS](https://emboss.sourceforge.net/docs/)                | `em_cons`                   | 6.6.0.0 |
| [Primer3](https://primer3.org/manual.html#invokingPrimer3)    | `primer3_core`              | 1.1.4   |


Python version employed: 3.11.13.  
Proyect python dependencies are listed in [**requirements.txt** ](./requirements.txt) file.  
Biopython and Snakemake drives most of the pipeline logic.

## Configuration:

The [**config/config.yml**](./config/config.yml) file contains parameters for the [*Viola*](https://en.wikipedia.org/wiki/Viola_(plant)) and [matK](https://en.wikipedia.org/wiki/Maturase_K) example.  
See [**notebooks/primer_evaluation.ipynb**](./notebooks/primer_evaluation.ipynb) for sample results and an evaluation walkthrough.

### config.yml fields
```yml
entrez:
  genes: # list[str] of gene names
  organisms: # list[str] of organism names
  min_len: # minimum sequence length (bp)
  max_len:  # maximum sequence length (bp)
primer3:
  PRIMER_OPT_SIZE:  # optimal primer length (bp)
  PRIMER_MIN_SIZE:  # minimum primer length (bp)
  PRIMER_MAX_SIZE:  # maximum primer length (bp)
  PRIMER_PRODUCT_SIZE_RANGE:  # range of amplicon size, e.g., "100-300"

```
### **.env** file and `email` environment variable
Include a .env file in the root project directory with the following content:
```dotenv
email="yourmail@example.com"
```
This email is required by NCBI when making Entrez queries.

## Data Pipeline
The pipeline designs candidate primers for a specific gene region in a given taxonomic group using NCBI resources and bioinformatics tools. It automates downloading sequences, multiple alignment, consensus generation, and Primer3 input preparation.
The workflow is defined in the Snakefile (Snakemake).

The DAG (Directed Acyclic Graph) of the pipeline is illustrated below:

<p align="center">
  <img src="grafo.png" alt="Pipeline DAG" width="120"><br>
  <small>Example for Viola and matK</small>
</p>

#### Search Logs
Information about each sequences retrieved for each search term (organism x gene) is saved as a separate `.log` file in the **logs/** directory.

### Running the Pipeline
From the root directory of the project, run:
```bash
snakemake --core 1
```
See ([official Snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html)) for details.


## Primer Evaluation 
### See: [notebooks/primer_evaluation.ipynb](./notebooks/primer_evaluation.ipynb)
The pipeline outputs a set of primers with scoring information.
This notebook demonstrates how to evaluate and interpret the results using Biopython and NCBI resources. It includes an example evaluating the best primer pair found for the *Viola* matK gene.